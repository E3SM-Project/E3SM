!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!> This module was derived from the original column physics driver by Adrian K. Turner and modified by the Icepack merge team.
!>
!
!-----------------------------------------------------------------------

module seaice_icepack

  use icepack_intfc, only: &
      icepack_warnings_aborted

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_timer
  use mpas_log, only: mpas_log_write, mpas_log_info

  use seaice_error

  implicit none

  private
  save

  public :: &
       seaice_init_icepack_physics_package_parameters, &
       seaice_init_icepack_physics_package_variables, &
       seaice_icepack_predynamics_time_integration, &
       seaice_icepack_dynamics_time_integration, &
       seaice_icepack_postdynamics_time_integration, &
       seaice_init_icepack_shortwave, &
       seaice_icepack_aggregate, &
       seaice_icepack_aggregate_simple, &
       seaice_icepack_initial_air_drag_coefficient, &
       seaice_icepack_reinitialize_fluxes, &
       seaice_icepack_reinitialize_diagnostics_thermodynamics, &
       seaice_icepack_reinitialize_diagnostics_dynamics, &
       seaice_icepack_reinitialize_diagnostics_bgc, &
       seaice_icepack_coupling_prep, &
       seaice_icepack_finalize, &
       seaice_icepack_init_trcr, &
       seaice_icepack_init_ocean_conc, &
       seaice_icepack_enthalpy_ice, &
       seaice_init_icepack_constants, &
       seaice_icepack_write_warnings

  ! tracer object
  type, private :: ciceTracerObjectType

     !-----------------------------------------------------------------------
     ! base tracer object
     !-----------------------------------------------------------------------

     ! length of tracer array
     integer :: nTracers   !ntrcr

     ! number of base tracers
     integer :: nBaseTracers = 3

     ! maximum number of ancestor tracers
     integer :: nMaxAncestorTracers = 2

     ! index of the parent tracer
     integer, dimension(:), allocatable :: parentIndex ! trcr_depend

     ! first ancestor type mask
     real(kind=RKIND), dimension(:,:), allocatable :: firstAncestorMask !trcr_base

     ! indices of ancestor tracers excluding base tracer
     integer, dimension(:,:), allocatable :: ancestorIndices ! nt_strata

     ! number of ancestor tracers excluding base tracer
     integer, dimension(:), allocatable :: ancestorNumber ! n_trcr_strata

     !-----------------------------------------------------------------------
     ! physics
     !-----------------------------------------------------------------------

     ! indexes of physics tracers in tracer array
     integer :: &
          index_surfaceTemperature, &    ! nt_Tsfc
          index_iceEnthalpy, &           ! nt_qice
          index_snowEnthalpy, &          ! nt_qsno
          index_iceSalinity, &           ! nt_sice
          index_iceAge, &                ! nt_iage
          index_firstYearIceArea, &      ! nt_FY
          index_levelIceArea, &          ! nt_alvl
          index_levelIceVolume, &        ! nt_vlvl
          index_pondArea, &              ! nt_apnd
          index_pondDepth, &             ! nt_hpnd
          index_pondLidThickness, &      ! nt_ipnd
          index_aerosols, &              ! nt_aero
          index_snowIceMass, &           ! nt_smice
          index_snowLiquidMass, &        ! nt_smliq
          index_snowGrainRadius, &       ! nt_rsnw
          index_snowDensity              ! nt_rhos

     !-----------------------------------------------------------------------
     ! biogeochemistry
     !-----------------------------------------------------------------------

     ! length of tracer array not including biology and biology related tracers
     integer :: nTracersNotBio   !ntrcr - ntrcr(bio)

     ! length of bio tracer array (does not include biology related tracers)
     integer :: nBioTracers   !nbtrcr

     ! number of bio tracers used (does not include brine or mobilefraction)
     integer :: nBioTracersLayer   !nltrcr

     ! length of shortwave bio tracer array (for aerosols and chlorophyll)
     integer :: nBioTracersShortwave   !nbtrcr_sw

     ! length of bio indices
     integer :: &
          nAlgaeIndex, &
          nAlgalCarbonIndex, &
          nAlgalChlorophyllIndex, &
          nDOCIndex, &
          nDONIndex, &
          nDICIndex, &
          nDissolvedIronIndex, &
          nParticulateIronIndex, &
          nzAerosolsIndex

     ! indexes of BGC tracers in tracer array
     integer :: &
          index_brineFraction, &         ! nt_fbri
          index_nitrateConc, &           ! nt_bgc_Nit
          index_ammoniumConc, &          ! nt_bgc_Am
          index_silicateConc, &          ! nt_bgc_Sil
          index_DMSPpConc, &             ! nt_bgc_DMSPp
          index_DMSPdConc, &             ! nt_bgc_DMSPd
          index_DMSConc, &               ! nt_bgc_DMS
          index_nonreactiveConc, &       ! nt_bgc_PON
          index_humicsConc, &            ! nt_bgc_hum
          index_mobileFraction, &        ! nt_zbgc_frac
          index_verticalSalinity, &      ! nt_bgc_S
          index_chlorophyllShortwave, &  ! nlt_chl_sw
          index_nitrateConcLayer, &      ! nlt_bgc_Nit
          index_ammoniumConcLayer, &     ! nlt_bgc_Am
          index_silicateConcLayer, &     ! nlt_bgc_Sil
          index_DMSPpConcLayer, &        ! nlt_bgc_DMSPp
          index_DMSPdConcLayer, &        ! nlt_bgc_DMSPd
          index_DMSConcLayer, &          ! nlt_bgc_DMS
          index_nonreactiveConcLayer, &  ! nlt_bgc_PON
          index_humicsConcLayer          ! nlt_bgc_hum

     ! indices of bio tracers with types in tracer array
     integer, dimension(:), allocatable :: &
          index_algaeConc, &                       ! nt_bgc_N
          index_algalCarbon, &                     ! nt_bgc_C
          index_algalChlorophyll, &                ! nt_bgc_chl
          index_DOCConc, &                         ! nt_bgc_DOC
          index_DONConc, &                         ! nt_bgc_DON
          index_DICConc, &                         ! nt_bgc_DIC
          index_dissolvedIronConc, &               ! nt_bgc_Fed
          index_particulateIronConc, &             ! nt_bgc_Fep
          index_verticalAerosolsConc, &            ! nt_zaero
          index_algaeConcLayer, &                  ! nlt_bgc_N
          index_algalCarbonLayer, &                ! nlt_bgc_C
          index_algalChlorophyllLayer, &           ! nlt_bgc_chl
          index_DOCConcLayer, &                    ! nlt_bgc_DOC
          index_DONConcLayer, &                    ! nlt_bgc_DON
          index_DICConcLayer, &                    ! nlt_bgc_DIC
          index_dissolvedIronConcLayer, &          ! nlt_bgc_Fed
          index_particulateIronConcLayer, &        ! nlt_bgc_Fep
          index_verticalAerosolsConcLayer, &       ! nlt_zaero
          index_verticalAerosolsConcShortwave, &   ! nlt_zaero_sw
          index_LayerIndexToDataArray, &           ! relates nlt to data array
          index_LayerIndexToBioIndex               ! relates nlt to nt

  end type ciceTracerObjectType

  type(ciceTracerObjectType), private :: ciceTracerObject

  real(kind=RKIND), dimension(:,:), allocatable :: &
       tracerArrayCategory
!$omp threadprivate(tracerArrayCategory)

  real(kind=RKIND), dimension(:), allocatable :: &
       tracerArrayCell

  integer, parameter :: strKINDWarnings = strKIND

contains

!-----------------------------------------------------------------------
! Initialize Column Physics Package
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_icepack_physics_package_parameters
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_icepack_physics_package_parameters(domain)

    type(domain_type), intent(inout) :: domain

    logical, pointer :: &
         config_use_column_physics

    call MPAS_pool_get_config(domain % configs, "config_use_column_physics", config_use_column_physics)
    if (config_use_column_physics) then

       ! set non activated variable pointers to other memory
       call init_icepack_non_activated_pointers(domain)

       ! initialize the column package tracer object
       call init_column_tracer_object(domain, ciceTracerObject)

       ! initialize the column package parameters
       call init_icepack_package_parameters(domain, ciceTracerObject)

    endif

  end subroutine seaice_init_icepack_physics_package_parameters

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_icepack_physics_package_variables
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_icepack_physics_package_variables(domain, clock)

    type(domain_type), intent(inout) :: domain

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    logical, pointer :: &
         config_use_column_physics, &
         config_do_restart, &
         config_use_column_biogeochemistry, &
         config_use_column_shortwave, &
         config_use_column_snow_tracers, &
         config_use_zaerosols

    call MPAS_pool_get_config(domain % configs, "config_use_column_physics", config_use_column_physics)
    call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_zaerosols", config_use_zaerosols)

    if (config_use_column_physics) then

       ! initialize level ice tracers
       call init_column_level_ice_tracers(domain)

       ! initialize the itd thickness classes
       call init_column_itd(domain)

       ! initialize thermodynamic tracer profiles
       call init_column_thermodynamic_profiles(domain)

       ! initialize biogoechemistry profiles
       if (config_use_column_biogeochemistry .or. config_use_zaerosols) &
            call init_column_biogeochemistry_profiles(domain, ciceTracerObject)

       ! history variables
       call init_column_history_variables(domain)

       ! snow
       call MPAS_pool_get_config(domain % configs, "config_use_column_snow_tracers", config_use_column_snow_tracers)
       if (config_use_column_snow_tracers) &
            call init_icepack_snow_tracers(domain)

       ! shortwave
       call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)
       call MPAS_pool_get_config(domain % configs, "config_use_column_shortwave", config_use_column_shortwave)
       if (config_do_restart .and. config_use_column_shortwave) &
            call seaice_init_icepack_shortwave(domain, clock)

    endif

  end subroutine seaice_init_icepack_physics_package_variables

!-----------------------------------------------------------------------
! column package initialization routine wrappers
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_itd
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_itd(domain)

    use icepack_intfc, only: &
         icepack_init_itd

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         initial

    real(kind=RKIND), dimension(:), pointer :: &
         categoryThicknessLimits

    integer, pointer :: &
         nCategories

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nCategories", nCategories)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "initial", initial)

       call MPAS_pool_get_array(initial, "categoryThicknessLimits", categoryThicknessLimits)

       call icepack_init_itd(&
            categoryThicknessLimits)

       call seaice_icepack_write_warnings(icepack_warnings_aborted())

       block => block % next
    end do

  end subroutine init_column_itd

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_thermodynamic_profiles
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_thermodynamic_profiles(domain)

    use icepack_intfc, only: &
         icepack_init_salinity, &
         icepack_liquidus_temperature

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         initial

    integer, pointer :: &
         nCellsSolve, &
         nIceLayers

    integer :: &
         iCell, &
         iIceLayer

    real(kind=RKIND), dimension(:), allocatable :: &
         initialSalinityProfileVertical

    real(kind=RKIND), dimension(:,:), pointer :: &
         initialSalinityProfile, &
         initialMeltingTemperatureProfile

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nIceLayers", nIceLayers)

       allocate(initialSalinityProfileVertical(1:nIceLayers+1))

       call icepack_init_salinity(&
            initialSalinityProfileVertical)
       call seaice_icepack_write_warnings(icepack_warnings_aborted())

       call MPAS_pool_get_subpool(block % structs, "initial", initial)

       call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)
       call MPAS_pool_get_array(initial, "initialMeltingTemperatureProfile", initialMeltingTemperatureProfile)

       do iCell = 1, nCellsSolve
          do iIceLayer = 1, nIceLayers + 1

             ! these profiles are not used by mushy
             initialSalinityProfile(iIceLayer,iCell)           = initialSalinityProfileVertical(iIceLayer)
             initialMeltingTemperatureProfile(iIceLayer,iCell) = &
                  icepack_liquidus_temperature(initialSalinityProfileVertical(iIceLayer))

          enddo ! iIceLayer
       enddo ! iCell

       deallocate(initialSalinityProfileVertical)

       block => block % next
    end do

  end subroutine init_column_thermodynamic_profiles

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_icepack_snow_tracers
!
!> \brief Initializes snow physics tracers
!>
!> \author Nicole Jeffery, LANL
!> \date 1 April 2017; modified for icepack July 2023
!> \details
!>
!> The following snow tracers are initialized:
!> 1) Snow liquid content (used to compute wet metamorphism of snow grain, modifies
!> liquid content of ponds, used in calculation of effective snow density due to content)
!> 2) Snow ice content (Used in calculation of effective snow density due to content)
!> 3) Effective snow density (both content and compaction are included. May be used for snow grain aging)
!> 4) Snow grain radius (used in radiative transfer calculations)
!
!-----------------------------------------------------------------------

  subroutine init_icepack_snow_tracers(domain)

    use seaice_constants, only: &
         seaicePuny, &
         seaiceDensitySnow

    use icepack_intfc, only: &
         icepack_init_snow

    type(domain_type), intent(in) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate, &
         snow

    logical, pointer :: &
         config_use_effective_snow_density, &
         config_use_snow_grain_radius, &
         config_do_restart_snow_density, &
         config_do_restart_snow_grain_radius

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         snowIceMass, &
         snowLiquidMass, &
         snowDensity, &
         snowVolumeCategory, &
         snowGrainRadius

    real(kind=RKIND), dimension(:,:), pointer :: &
         snowMeltMassCategory

    real(kind=RKIND), dimension(:), pointer :: &
         snowDensityViaContent, &
         snowDensityViaCompaction, &
         snowMeltMassCell, &
         snowVolumeCell

    real(kind=RKIND), pointer :: &
         config_fallen_snow_radius, &
         config_new_snow_density

    integer, pointer :: &
         nCellsSolve, &
         nSnowLayers, &
         nCategories

    integer :: &
         iCell, &
         iSnowLayer, &
         iCategory

    call icepack_init_snow()

    call MPAS_pool_get_config(domain % configs, "config_use_effective_snow_density", config_use_effective_snow_density)
    call MPAS_pool_get_config(domain % configs, "config_use_snow_grain_radius", config_use_snow_grain_radius)
    call MPAS_pool_get_config(domain % configs, "config_fallen_snow_radius", config_fallen_snow_radius)
    call MPAS_pool_get_config(domain % configs, "config_new_snow_density", config_new_snow_density)
    call MPAS_pool_get_config(domain % configs, "config_do_restart_snow_density", config_do_restart_snow_density)
    call MPAS_pool_get_config(domain % configs, "config_do_restart_snow_grain_radius", config_do_restart_snow_grain_radius)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
       call MPAS_pool_get_subpool(block % structs, "snow", snow)

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nCategories", nCategories)
       call MPAS_pool_get_dimension(block % dimensions, "nSnowLayers", nSnowLayers)

       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
       call MPAS_pool_get_array(tracers_aggregate, "snowVolumeCell", snowVolumeCell)
       call MPAS_pool_get_array(snow, "snowMeltMassCategory", snowMeltMassCategory)
       call MPAS_pool_get_array(snow, "snowMeltMassCell", snowMeltMassCell)
       call MPAS_pool_get_array(snow, "snowDensityViaContent", snowDensityViaContent)

       snowMeltMassCategory(:,:)   = 0.0_RKIND
       snowMeltMassCell(:)         = 0.0_RKIND
       snowDensityViaContent(:)    = seaiceDensitySnow

       if (config_use_effective_snow_density) then

          call MPAS_pool_get_array(snow, "snowDensityViaCompaction", snowDensityViaCompaction)
          call MPAS_pool_get_array(tracers, "snowDensity", snowDensity, 1)

          if (.not. config_do_restart_snow_density) then
             snowDensity(:,:,:)          = config_new_snow_density
             do iCell = 1, nCellsSolve
                if (snowVolumeCell(iCell) .gt. seaicePuny) then
                   snowDensityViaCompaction(iCell) = config_new_snow_density
                else
                   snowDensityViaCompaction(iCell) = 0.0_RKIND
                endif
             enddo
          else
             snowDensityViaCompaction(:) = 0.0_RKIND

             do iCell = 1, nCellsSolve

                do iCategory = 1, nCategories
                   if (snowVolumeCategory(1,iCategory,iCell) .gt. 0.0_RKIND) then
                      do iSnowLayer = 1, nSnowLayers
                         snowDensityViaCompaction(iCell) = snowDensityViaCompaction(iCell) &
                            + snowVolumeCategory(1,iCategory,iCell) * &
                            snowDensity(iSnowLayer,iCategory,iCell)
                      enddo !iSnowLayer
                   endif !snowVolumeCategory
                enddo !iCategory
                if (snowVolumeCell(iCell) .gt.  seaicePuny) then
                   snowDensityViaCompaction(iCell) = snowDensityViaCompaction(iCell)/ &
                      (snowVolumeCell(iCell) * real(nSnowLayers,kind=RKIND))
                else
                   snowDensityViaCompaction(iCell) = 0.0_RKIND
                endif !snowVolumeCell

             enddo ! iCell
          endif ! config_do_restart_snow_density
       endif !config_use_effective_snow_density

       if (config_use_snow_grain_radius) then

          call MPAS_pool_get_array(tracers, "snowIceMass", snowIceMass, 1)
          call MPAS_pool_get_array(tracers, "snowLiquidMass", snowLiquidMass, 1)
          call MPAS_pool_get_array(tracers, "snowGrainRadius", snowGrainRadius, 1)

          if (.not. config_do_restart_snow_grain_radius) then

             snowGrainRadius(:,:,:)      = config_fallen_snow_radius
             snowIceMass(:,:,:)          = seaiceDensitySnow
             snowLiquidMass(:,:,:)       = 0.0_RKIND
             snowDensityViaContent(:)    = seaiceDensitySnow

          else
             do iCell = 1, nCellsSolve
                snowDensityViaContent(iCell) = 0.0_RKIND
                do iCategory = 1, nCategories
                   if (snowVolumeCategory(1,iCategory,iCell) .gt. 0.0_RKIND) then
                      do iSnowLayer = 1, nSnowLayers
                         snowDensityViaContent(iCell) = snowDensityViaContent(iCell) &
                             + snowVolumeCategory(1,iCategory,iCell) * &
                             (snowIceMass(iSnowLayer,iCategory,iCell) + &
                             snowLiquidMass(iSnowLayer,iCategory,iCell))
                       enddo !iSnowLayer
                   endif !snowVolumeCategory
                enddo !iCategory
                if (snowVolumeCell(iCell) .gt.  seaicePuny) then
                   snowDensityViaContent(iCell) = snowDensityViaContent(iCell)/ &
                       (snowVolumeCell(iCell) * real(nSnowLayers,kind=RKIND))  !!!CHECK THIS!!!
                else
                   snowDensityViaContent(iCell)    = seaiceDensitySnow
                endif !snowVolumeCell
             enddo ! iCell
          endif ! config_do_restart_snow_grain_radius
        endif ! config_use_snow_grain_radius
        block => block % next
     end do

  end subroutine init_icepack_snow_tracers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_icepack_shortwave
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_icepack_shortwave(domain, clock)

    use icepack_intfc, only: &
         icepack_init_orbit, &
         icepack_init_radiation

    use seaice_constants, only: &
         seaicePuny

    type(domain_type), intent(inout) :: domain

    type(MPAS_clock_type), intent(in) :: clock

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         shortwave, &
         atmos_coupling

    real(kind=RKIND), dimension(:), pointer :: &
         solarZenithAngleCosine, &
         albedoVisibleDirectCell, &
         albedoVisibleDiffuseCell, &
         albedoIRDirectCell, &
         albedoIRDiffuseCell, &
         albedoVisibleDirectArea, &
         albedoVisibleDiffuseArea, &
         albedoIRDirectArea, &
         albedoIRDiffuseArea, &
         bareIceAlbedoCell, &
         snowAlbedoCell, &
         pondAlbedoCell, &
         effectivePondAreaCell, &
         shortwaveScalingFactor, &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown

    real(kind=RKIND), dimension(:,:), pointer :: &
         albedoVisibleDirectCategory, &
         albedoVisibleDiffuseCategory, &
         albedoIRDirectCategory, &
         albedoIRDiffuseCategory, &
         bareIceAlbedoCategory, &
         snowAlbedoCategory, &
         pondAlbedoCategory, &
         effectivePondAreaCategory

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    integer :: &
         iCell, &
         iCategory

    character(len=strKIND), pointer :: &
         config_shortwave_type

    logical, pointer :: &
         config_do_restart

    call icepack_init_radiation()
    call seaice_icepack_write_warnings(icepack_warnings_aborted())

    call MPAS_pool_get_config(domain % configs, "config_shortwave_type", config_shortwave_type)
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    if (trim(config_shortwave_type(1:4)) == "dEdd") then

       call icepack_init_orbit()
       call seaice_icepack_write_warnings(icepack_warnings_aborted())

    endif

    call column_radiation(&
         domain, &
         clock, &
         .true.)

    ! other shortwave initialization
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)

       call MPAS_pool_get_array(shortwave, "solarZenithAngleCosine", solarZenithAngleCosine)

       call MPAS_pool_get_array(shortwave, "albedoVisibleDirectCell", albedoVisibleDirectCell)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDiffuseCell", albedoVisibleDiffuseCell)
       call MPAS_pool_get_array(shortwave, "albedoIRDirectCell", albedoIRDirectCell)
       call MPAS_pool_get_array(shortwave, "albedoIRDiffuseCell", albedoIRDiffuseCell)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDirectCategory", albedoVisibleDirectCategory)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDiffuseCategory", albedoVisibleDiffuseCategory)
       call MPAS_pool_get_array(shortwave, "albedoIRDirectCategory", albedoIRDirectCategory)
       call MPAS_pool_get_array(shortwave, "albedoIRDiffuseCategory", albedoIRDiffuseCategory)

       call MPAS_pool_get_array(shortwave, "albedoVisibleDirectArea", albedoVisibleDirectArea)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDiffuseArea", albedoVisibleDiffuseArea)
       call MPAS_pool_get_array(shortwave, "albedoIRDirectArea", albedoIRDirectArea)
       call MPAS_pool_get_array(shortwave, "albedoIRDiffuseArea", albedoIRDiffuseArea)

       call MPAS_pool_get_array(shortwave, "bareIceAlbedoCell", bareIceAlbedoCell)
       call MPAS_pool_get_array(shortwave, "snowAlbedoCell", snowAlbedoCell)
       call MPAS_pool_get_array(shortwave, "pondAlbedoCell", pondAlbedoCell)

       call MPAS_pool_get_array(shortwave, "bareIceAlbedoCategory", bareIceAlbedoCategory)
       call MPAS_pool_get_array(shortwave, "snowAlbedoCategory", snowAlbedoCategory)
       call MPAS_pool_get_array(shortwave, "pondAlbedoCategory", pondAlbedoCategory)

       call MPAS_pool_get_array(shortwave, "effectivePondAreaCell", effectivePondAreaCell)
       call MPAS_pool_get_array(shortwave, "effectivePondAreaCategory", effectivePondAreaCategory)

       call MPAS_pool_get_array(shortwave, "shortwaveScalingFactor", shortwaveScalingFactor)

       call MPAS_pool_get_array(atmos_coupling, "shortwaveVisibleDirectDown", shortwaveVisibleDirectDown)
       call MPAS_pool_get_array(atmos_coupling, "shortwaveVisibleDiffuseDown", shortwaveVisibleDiffuseDown)
       call MPAS_pool_get_array(atmos_coupling, "shortwaveIRDirectDown", shortwaveIRDirectDown)
       call MPAS_pool_get_array(atmos_coupling, "shortwaveIRDiffuseDown", shortwaveIRDiffuseDown)

       do iCell = 1, nCellsSolve

          albedoVisibleDirectCell(iCell)  = 0.0_RKIND
          albedoVisibleDiffuseCell(iCell) = 0.0_RKIND
          albedoIRDirectCell(iCell)       = 0.0_RKIND
          albedoIRDiffuseCell(iCell)      = 0.0_RKIND

          do iCategory = 1, nCategories

             ! aggregate albedos
             if (iceAreaCategory(1,iCategory,iCell) > seaicePuny) then

                albedoVisibleDirectCell(iCell)  = albedoVisibleDirectCell(iCell)  + &
                     albedoVisibleDirectCategory(iCategory,iCell)  * iceAreaCategory(1,iCategory,iCell)
                albedoVisibleDiffuseCell(iCell) = albedoVisibleDiffuseCell(iCell) + &
                     albedoVisibleDiffuseCategory(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)
                albedoIRDirectCell(iCell)       = albedoIRDirectCell(iCell)       + &
                     albedoIRDirectCategory(iCategory,iCell)       * iceAreaCategory(1,iCategory,iCell)
                albedoIRDiffuseCell(iCell)      = albedoIRDiffuseCell(iCell)      + &
                     albedoIRDiffuseCategory(iCategory,iCell)      * iceAreaCategory(1,iCategory,iCell)

                if (solarZenithAngleCosine(iCell) > seaicePuny) then ! sun above horizon

                   bareIceAlbedoCell(iCell) = bareIceAlbedoCell(iCell) + &
                        bareIceAlbedoCategory(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)
                   snowAlbedoCell(iCell)    = snowAlbedoCell(iCell)    + &
                        snowAlbedoCategory(iCategory,iCell)    * iceAreaCategory(1,iCategory,iCell)
                   pondAlbedoCell(iCell)    = pondAlbedoCell(iCell)    + &
                        pondAlbedoCategory(iCategory,iCell)    * iceAreaCategory(1,iCategory,iCell)

                endif

                effectivePondAreaCell(iCell) = effectivePondAreaCell(iCell) + &
                     effectivePondAreaCategory(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)

             endif

          enddo ! iCategory

          ! Store grid box mean albedos and fluxes before scaling by aice
          albedoVisibleDirectArea(iCell)  = albedoVisibleDirectCell(iCell)
          albedoVisibleDiffuseArea(iCell) = albedoVisibleDiffuseCell(iCell)
          albedoIRDirectArea(iCell)       = albedoIRDirectCell(iCell)
          albedoIRDiffuseArea(iCell)      = albedoIRDiffuseCell(iCell)

          ! Save net shortwave for scaling factor in scale_factor
          if (.not. config_do_restart) then
             shortwaveScalingFactor(iCell) = &
                  shortwaveVisibleDirectDown(iCell)  * (1.0_RKIND - albedoVisibleDirectArea(iCell)) + &
                  shortwaveVisibleDiffuseDown(iCell) * (1.0_RKIND - albedoVisibleDiffuseArea(iCell)) + &
                  shortwaveIRDirectDown(iCell)       * (1.0_RKIND - albedoIRDirectArea(iCell)) + &
                  shortwaveIRDiffuseDown(iCell)      * (1.0_RKIND - albedoIRDiffuseArea(iCell))
          endif

       enddo ! iCell

       block => block % next
    end do

  end subroutine seaice_init_icepack_shortwave

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_thermodynamic_tracers
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_thermodynamic_tracers(domain)

    use icepack_intfc, only: &
         icepack_init_enthalpy

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         atmos_coupling, &
         ocean_coupling, &
         initial

    real(kind=RKIND), dimension(:), pointer :: &
         airTemperature, &
         seaFreezingTemperature

    real(kind=RKIND), dimension(:,:), pointer :: &
         initialSalinityProfile, &
         initialMeltingTemperatureProfile

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         surfaceTemperature, &
         iceEnthalpy, &
         snowEnthalpy

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    integer :: &
         iCell, &
         iCategory

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
       call MPAS_pool_get_subpool(block % structs, "initial", initial)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

       call MPAS_pool_get_array(atmos_coupling, "airTemperature", airTemperature)

       call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)

       call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)
       call MPAS_pool_get_array(initial, "initialMeltingTemperatureProfile", initialMeltingTemperatureProfile)

       call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
       call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy, 1)
       call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy, 1)

       do iCell = 1, nCellsSolve
          do iCategory = 1, nCategories

             call icepack_init_enthalpy(&
                  airTemperature(iCell), &
                  seaFreezingTemperature(iCell), &
                  initialSalinityProfile(:,iCell), &
                  initialMeltingTemperatureProfile(:,iCell), &
                  surfaceTemperature(1,iCategory,iCell), &
                  iceEnthalpy(:,iCategory,iCell), &
                  snowEnthalpy(:,iCategory,iCell))

          enddo ! iCategory
       enddo ! iCell

       block => block % next
    end do

    call seaice_icepack_write_warnings(icepack_warnings_aborted())

  end subroutine init_column_thermodynamic_tracers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_level_ice_tracers
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_level_ice_tracers(domain)

    type(domain_type), intent(inout) :: domain

    logical, pointer :: &
         config_use_level_ice, &
         config_do_restart

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         tracers

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         levelIceArea, &
         levelIceVolume

    call MPAS_pool_get_config(domain % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    if (config_use_level_ice .and. .not. config_do_restart) then

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

          call MPAS_pool_get_array(tracers, "levelIceArea", levelIceArea, 1)
          call MPAS_pool_get_array(tracers, "levelIceVolume", levelIceVolume, 1)

          levelIceArea = 1.0_RKIND
          levelIceVolume = 1.0_RKIND

          block => block % next
       end do

    endif

  end subroutine init_column_level_ice_tracers

!-----------------------------------------------------------------------
! finalize
!-----------------------------------------------------------------------
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack_finalize
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_icepack_finalize(domain)

    type(domain_type), intent(inout) :: domain

    call finalize_icepack_non_activated_pointers(domain)

  end subroutine seaice_icepack_finalize

!-----------------------------------------------------------------------
! runtime
!-----------------------------------------------------------------------
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack_predynamics_time_integration
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_icepack_predynamics_time_integration(domain, clock)

    type(domain_type), intent(inout) :: domain

    type(MPAS_clock_type), intent(in) :: clock

    logical, pointer :: &
         config_use_column_physics, &
         config_use_column_shortwave, &
         config_use_column_vertical_thermodynamics, &
         config_use_column_biogeochemistry, &
         config_use_column_itd_thermodynamics, &
         config_calc_surface_temperature, &
         config_use_vertical_tracers, &
         config_use_zaerosols

    real(kind=RKIND), pointer :: &
         config_dt

    call MPAS_pool_get_config(domain % configs, "config_use_column_physics", config_use_column_physics)

    if (config_use_column_physics) then

       call MPAS_pool_get_config(domain % configs, "config_use_column_shortwave", config_use_column_shortwave)
       call MPAS_pool_get_config(domain % configs, "config_use_column_vertical_thermodynamics", &
                                                    config_use_column_vertical_thermodynamics)
       call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
       call MPAS_pool_get_config(domain % configs, "config_use_zaerosols", config_use_zaerosols)
       call MPAS_pool_get_config(domain % configs, "config_use_column_itd_thermodynamics", config_use_column_itd_thermodynamics)
       call MPAS_pool_get_config(domain % configs, "config_calc_surface_temperature", config_calc_surface_temperature)
       call MPAS_pool_get_config(domain % configs, "config_use_vertical_tracers", config_use_vertical_tracers)

       call MPAS_pool_get_config(domain % configs, "config_dt", config_dt)

       !-----------------------------------------------------------------
       ! Scale radiation fields
       !-----------------------------------------------------------------

       call mpas_timer_start("Column shortwave prep")
       if (config_use_column_shortwave .and. config_calc_surface_temperature) &
            call column_prep_radiation(domain)
       call mpas_timer_stop("Column shortwave prep")

       !-----------------------------------------------------------------
       ! Vertical thermodynamics
       !-----------------------------------------------------------------

       call mpas_timer_start("Column vertical thermodynamics")
       if (config_use_column_vertical_thermodynamics) &
            call column_vertical_thermodynamics(domain, clock)
       call mpas_timer_stop("Column vertical thermodynamics")

       !-----------------------------------------------------------------
       ! Biogeochemistry
       !-----------------------------------------------------------------

       call mpas_timer_start("Column biogeochemistry")
       if (config_use_column_biogeochemistry .or. config_use_zaerosols) &
            call column_biogeochemistry(domain)
       call mpas_timer_stop("Column biogeochemistry")

       !-----------------------------------------------------------------
       ! ITD thermodynamics
       !-----------------------------------------------------------------

       call mpas_timer_start("Column ITD thermodynamics")
       if (config_use_column_itd_thermodynamics) &
            call column_itd_thermodynamics(domain, clock)
       call mpas_timer_stop("Column ITD thermodynamics")

       !-----------------------------------------------------------------
       ! Update the aggregated state variables
       !-----------------------------------------------------------------

       call mpas_timer_start("Column predyn update state")
       call seaice_column_update_state(domain, "thermodynamics", config_dt, config_dt)
       call mpas_timer_stop("Column predyn update state")

       !-----------------------------------------------------------------
       ! Separate vertical snow and ice tracers for advection
       !-----------------------------------------------------------------

       call mpas_timer_start("Column separate snow/ice tracers")
       if (config_use_vertical_tracers) &
            call column_separate_snow_ice_tracers(domain)
       call mpas_timer_stop("Column separate snow/ice tracers")

    endif

  end subroutine seaice_icepack_predynamics_time_integration

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack_dynamics_time_integration
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_icepack_dynamics_time_integration(domain, clock)

    type(domain_type), intent(inout) :: domain

    type(MPAS_clock_type), intent(in) :: clock

    logical, pointer :: &
         config_use_column_physics, &
         config_use_column_ridging, &
         config_use_vertical_tracers

    type(MPAS_pool_type), pointer :: &
         velocitySolver

    real(kind=RKIND), pointer :: &
         dynamicsTimeStep

    call MPAS_pool_get_config(domain % configs, "config_use_column_physics", config_use_column_physics)

    if (config_use_column_physics) then

       call MPAS_pool_get_config(domain % configs, "config_use_column_ridging", config_use_column_ridging)
       call MPAS_pool_get_config(domain % configs, "config_use_vertical_tracers", config_use_vertical_tracers)

       call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolver)
       call MPAS_pool_get_array(velocitySolver, "dynamicsTimeStep", dynamicsTimeStep)

       !-----------------------------------------------------------------
       ! Combine vertical snow and ice tracers
       !-----------------------------------------------------------------

       call mpas_timer_start("Column combine snow/ice tracers")
       if (config_use_vertical_tracers) &
            call column_combine_snow_ice_tracers(domain)
       call mpas_timer_stop("Column combine snow/ice tracers")

       !-----------------------------------------------------------------
       ! Ridging
       !-----------------------------------------------------------------

       call mpas_timer_start("Column ridging")
       if (config_use_column_ridging) &
            call column_ridging(domain)
       call mpas_timer_stop("Column ridging")

       !-----------------------------------------------------------------
       ! Update the aggregated state variables
       !-----------------------------------------------------------------

       call mpas_timer_start("Column update state")
       call seaice_column_update_state(domain, "transport", dynamicsTimeStep, 0.0_RKIND)
       call mpas_timer_stop("Column update state")

    else

       call mpas_timer_start("Column aggregate")
       call seaice_icepack_aggregate_simple(domain)
       call mpas_timer_stop("Column aggregate")

    endif

  end subroutine seaice_icepack_dynamics_time_integration

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack_postdynamics_time_integration
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_icepack_postdynamics_time_integration(domain, clock)

    type(domain_type), intent(inout) :: domain

    type(MPAS_clock_type), intent(in) :: clock

    logical, pointer :: &
         config_use_column_physics, &
         config_use_column_shortwave, &
         config_use_column_snow_tracers

    type(block_type), pointer :: &
         block

    call MPAS_pool_get_config(domain % configs, "config_use_column_physics", config_use_column_physics)

    if (config_use_column_physics) then

       call MPAS_pool_get_config(domain % configs, "config_use_column_shortwave", config_use_column_shortwave)
       call MPAS_pool_get_config(domain % configs, "config_use_column_snow_tracers", config_use_column_snow_tracers)

       !-----------------------------------------------------------------
       ! snow
       !-----------------------------------------------------------------

       call mpas_timer_start("Column snow")
       if (config_use_column_snow_tracers) &
            call column_snow(domain)
       call mpas_timer_stop("Column snow")

       !-----------------------------------------------------------------
       ! Shortwave radiation
       !-----------------------------------------------------------------

       call mpas_timer_start("Column shortwave")
       if (config_use_column_shortwave) &
            call column_radiation(domain, clock, .false.)
       call mpas_timer_stop("Column shortwave")

       !-----------------------------------------------------------------
       ! Coupling prep
       !-----------------------------------------------------------------

       call mpas_timer_start("Column coupling prep")
       call seaice_icepack_coupling_prep(domain)
       call mpas_timer_stop("Column coupling prep")

    endif

  end subroutine seaice_icepack_postdynamics_time_integration

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_vertical_thermodynamics
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine column_vertical_thermodynamics(domain, clock)

    use icepack_intfc, only: &
         icepack_step_therm1, &
         icepack_warnings_clear

    use seaice_constants, only: &
         seaicePuny

    use seaice_mesh, only: &
         seaice_interpolate_vertex_to_cell

    type(domain_type), intent(inout) :: domain

    type(MPAS_clock_type), intent(in) :: clock

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         icestate, &
         tracers, &
         tracers_aggregate, &
         velocity_solver, &
         atmos_coupling, &
         atmos_forcing, &
         alternative_atmos_forcing, &
         ocean_coupling, &
         drag, &
         melt_growth_rates, &
         atmos_fluxes, &
         ocean_fluxes, &
         shortwave, &
         ponds, &
         aerosols, &
         diagnostics, &
         snow, &
         boundary

    ! configs
    real(kind=RKIND), pointer :: &
         config_dt

    logical, pointer :: &
         config_use_aerosols, &
         config_use_prescribed_ice, &
         config_use_snow_liquid_ponds, &
         config_use_high_frequency_coupling

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
         uVelocityCell, &
         vvelocityCell, &
         uAirVelocity, &
         vAirVelocity, &
         windSpeed, &
         airLevelHeight, &
         airSpecificHumidity, &
         airDensity, &
         airTemperature, &
         atmosReferenceTemperature2m, &
         atmosReferenceHumidity2m, &
         atmosReferenceSpeed10m, &
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
         bottomTemperature, &
         oceanStressCellU, &
         oceanStressCellV, &
         freezingMeltingPotential, &
         lateralIceMeltRate, &
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
         evaporativeWaterFlux, &
         evaporativeWaterFluxSnow, &
         evaporativeWaterFluxIce, &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         oceanShortwaveFlux, &
         bottomConductiveFlux, &
         surfaceIceMelt, &
         basalIceMelt, &
         lateralIceMelt, &
         snowMelt, &
         congelation, &
         snowiceFormation, &
         snowThicknessChangeCell, &
         frazilFormation, &
         meltOnset, &
         freezeOnset, &
         snowIceInterfaceTemperature, &
         oceanHeatFluxIceBottom, &
         openWaterArea, &
         snowLossToLeads, &
         snowMeltMassCell

    real(kind=RKIND), dimension(:,:), pointer :: &
         iceAreaCategoryInitial, &
         iceVolumeCategoryInitial, &
         snowVolumeCategoryInitial, &
         surfaceShortwaveFlux, &
         interiorShortwaveFlux, &
         penetratingShortwaveFlux, &
         sensibleHeatFluxCategory, &
         latentHeatFluxCategory, &
         surfaceIceMeltCategory, &
         basalIceMeltCategory, &
         lateralIceMeltFractionCategory, &
         snowMeltCategory, &
         congelationCategory, &
         snowiceFormationCategory, &
         atmosAerosolFlux, &
         oceanAerosolFlux, &
         pondSnowDepthDifference, &
         pondLidMeltFluxFraction, &
         surfaceHeatFluxCategory, &
         surfaceConductiveFluxCategory, &
         bottomConductiveFluxCategory, &
         latentHeatFluxCouple, &
         sensibleHeatFluxCouple, &
         surfaceHeatFluxCouple, &
         surfaceConductiveFluxCouple, &
         snowThicknessChangeCategory, &
         snowMeltMassCategory, &
         snowRadiusInStandardRadiationSchemeCategory

    real(kind=RKIND), dimension(:,:,:), pointer :: &
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
         snowEnthalpy, &
         iceEnthalpy, &
         iceSalinity, &
         absorbedShortwaveIceLayer, &
         absorbedShortwaveSnowLayer, &
         snowScatteringAerosol, &
         snowBodyAerosol, &
         iceScatteringAerosol, &
         iceBodyAerosol, &
         snowIceMass, &
         snowLiquidMass, &
         snowGrainRadius

    integer, dimension(:), pointer :: &
         indexToCellID

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
         abortMessage, &
         abortLocation

    real(kind=RKIND) :: &
         dayOfYear

    ! day of year
    call get_day_of_year(clock, dayOfYear)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
       call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmos_forcing)
       call MPAS_pool_get_subpool(block % structs, "alternative_atmos_forcing", alternative_atmos_forcing)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
       call MPAS_pool_get_subpool(block % structs, "drag", drag)
       call MPAS_pool_get_subpool(block % structs, "melt_growth_rates", melt_growth_rates)
       call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmos_fluxes)
       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", ocean_fluxes)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "ponds", ponds)
       call MPAS_pool_get_subpool(block % structs, "aerosols", aerosols)
       call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnostics)
       call MPAS_pool_get_subpool(block % structs, "snow", snow)
       call MPAS_pool_get_subpool(block % structs, "boundary", boundary)

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)
       call MPAS_pool_get_config(block % configs, "config_use_prescribed_ice", config_use_prescribed_ice)
       call MPAS_pool_get_config(block % configs, "config_use_snow_liquid_ponds", config_use_snow_liquid_ponds)
       call MPAS_pool_get_config(block % configs, "config_use_high_frequency_coupling", config_use_high_frequency_coupling)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)
       call MPAS_pool_get_dimension(mesh, "nAerosols", nAerosols)

       call MPAS_pool_get_array(mesh, "latCell", latCell)
       call MPAS_pool_get_array(mesh, "indexToCellID", indexToCellID)

       call MPAS_pool_get_array(icestate, "iceAreaCellInitial", iceAreaCellInitial)
       call MPAS_pool_get_array(icestate, "iceAreaCategoryInitial", iceAreaCategoryInitial)
       call MPAS_pool_get_array(icestate, "iceVolumeCategoryInitial", iceVolumeCategoryInitial)
       call MPAS_pool_get_array(icestate, "snowVolumeCategoryInitial", snowVolumeCategoryInitial)
       call MPAS_pool_get_array(icestate, "openWaterArea", openWaterArea)

       call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(tracers_aggregate, "iceVolumeCell", iceVolumeCell)
       call MPAS_pool_get_array(tracers_aggregate, "snowVolumeCell", snowVolumeCell)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
       call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy, 1)
       call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy, 1)
       call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity, 1)
       call MPAS_pool_get_array(tracers, "levelIceArea", levelIceArea, 1)
       call MPAS_pool_get_array(tracers, "levelIceVolume", levelIceVolume, 1)
       call MPAS_pool_get_array(tracers, "pondArea", pondArea, 1)
       call MPAS_pool_get_array(tracers, "pondDepth", pondDepth, 1)
       call MPAS_pool_get_array(tracers, "pondLidThickness", pondLidThickness, 1)
       call MPAS_pool_get_array(tracers, "iceAge", iceAge, 1)
       call MPAS_pool_get_array(tracers, "firstYearIceArea", firstYearIceArea, 1)
       call MPAS_pool_get_array(tracers, "snowScatteringAerosol", snowScatteringAerosol, 1)
       call MPAS_pool_get_array(tracers, "snowBodyAerosol", snowBodyAerosol, 1)
       call MPAS_pool_get_array(tracers, "iceScatteringAerosol", iceScatteringAerosol, 1)
       call MPAS_pool_get_array(tracers, "iceBodyAerosol", iceBodyAerosol, 1)
       call MPAS_pool_get_array(tracers, "snowIceMass", snowIceMass, 1)
       call MPAS_pool_get_array(tracers, "snowLiquidMass", snowLiquidMass, 1)
       call MPAS_pool_get_array(tracers, "snowGrainRadius", snowGrainRadius, 1)

       call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocity_solver, "uVelocityCell", uVelocityCell)
       call MPAS_pool_get_array(velocity_solver, "vVelocityCell", vVelocityCell)
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
       call MPAS_pool_get_array(atmos_coupling, "atmosReferenceTemperature2m", atmosReferenceTemperature2m)
       call MPAS_pool_get_array(atmos_coupling, "atmosReferenceHumidity2m", atmosReferenceHumidity2m)
       call MPAS_pool_get_array(atmos_coupling, "atmosReferenceSpeed10m", atmosReferenceSpeed10m)

       call MPAS_pool_get_array(atmos_forcing, "windSpeed", windSpeed)

       call MPAS_pool_get_array(alternative_atmos_forcing, "latentHeatFluxCouple", latentHeatFluxCouple)
       call MPAS_pool_get_array(alternative_atmos_forcing, "sensibleHeatFluxCouple", sensibleHeatFluxCouple)
       call MPAS_pool_get_array(alternative_atmos_forcing, "surfaceHeatFluxCouple", surfaceHeatFluxCouple)
       call MPAS_pool_get_array(alternative_atmos_forcing, "surfaceConductiveFluxCouple", surfaceConductiveFluxCouple)
       call MPAS_pool_get_array(alternative_atmos_forcing, "airStressForcingU", airStressForcingU)
       call MPAS_pool_get_array(alternative_atmos_forcing, "airStressForcingV", airStressForcingV)

       call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTemperature", seaSurfaceTemperature)
       call MPAS_pool_get_array(ocean_coupling, "seaSurfaceSalinity", seaSurfaceSalinity)
       call MPAS_pool_get_array(ocean_coupling, "bottomTemperature", bottomTemperature)
       call MPAS_pool_get_array(ocean_coupling, "freezingMeltingPotential", freezingMeltingPotential)
       call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)

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

       call MPAS_pool_get_array(melt_growth_rates, "lateralIceMeltFractionCategory", lateralIceMeltFractionCategory)
       call MPAS_pool_get_array(melt_growth_rates, "lateralIceMeltRate", lateralIceMeltRate)
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
       call MPAS_pool_get_array(melt_growth_rates, "snowThicknessChangeCell", snowThicknessChangeCell)
       call MPAS_pool_get_array(melt_growth_rates, "snowThicknessChangeCategory", snowThicknessChangeCategory)
       call MPAS_pool_get_array(melt_growth_rates, "frazilFormation", frazilFormation)

       call MPAS_pool_get_array(atmos_fluxes, "surfaceHeatFlux", surfaceHeatFlux)
       call MPAS_pool_get_array(atmos_fluxes, "surfaceHeatFluxCategory", surfaceHeatFluxCategory)
       call MPAS_pool_get_array(atmos_fluxes, "surfaceConductiveFlux", surfaceConductiveFlux)
       call MPAS_pool_get_array(atmos_fluxes, "surfaceConductiveFluxCategory", surfaceConductiveFluxCategory)
       call MPAS_pool_get_array(atmos_fluxes, "longwaveUp", longwaveUp)
       call MPAS_pool_get_array(atmos_fluxes, "sensibleHeatFlux", sensibleHeatFlux)
       call MPAS_pool_get_array(atmos_fluxes, "sensibleHeatFluxCategory", sensibleHeatFluxCategory)
       call MPAS_pool_get_array(atmos_fluxes, "latentHeatFlux", latentHeatFlux)
       call MPAS_pool_get_array(atmos_fluxes, "latentHeatFluxCategory", latentHeatFluxCategory)
       call MPAS_pool_get_array(atmos_fluxes, "evaporativeWaterFlux", evaporativeWaterFlux)
       call MPAS_pool_get_array(atmos_fluxes, "evaporativeWaterFluxSnow", evaporativeWaterFluxSnow)
       call MPAS_pool_get_array(atmos_fluxes, "evaporativeWaterFluxIce", evaporativeWaterFluxIce)

       call MPAS_pool_get_array(ocean_fluxes, "oceanFreshWaterFlux", oceanFreshWaterFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanSaltFlux", oceanSaltFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanHeatFlux", oceanHeatFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanShortwaveFlux", oceanShortwaveFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanHeatFluxIceBottom", oceanHeatFluxIceBottom)
       call MPAS_pool_get_array(ocean_fluxes, "bottomConductiveFlux", bottomConductiveFlux)
       call MPAS_pool_get_array(ocean_fluxes, "bottomConductiveFluxCategory", bottomConductiveFluxCategory)

       call MPAS_pool_get_array(shortwave, "surfaceShortwaveFlux", surfaceShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "interiorShortwaveFlux", interiorShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "penetratingShortwaveFlux", penetratingShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "absorbedShortwaveFlux", absorbedShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "absorbedShortwaveIceLayer", absorbedShortwaveIceLayer)
       call MPAS_pool_get_array(shortwave, "absorbedShortwaveSnowLayer", absorbedShortwaveSnowLayer)
       call MPAS_pool_get_array(shortwave, "solarZenithAngleCosine", solarZenithAngleCosine)

       call MPAS_pool_get_array(aerosols, "atmosAerosolFlux", atmosAerosolFlux)
       call MPAS_pool_get_array(aerosols, "oceanAerosolFlux", oceanAerosolFlux)

       call MPAS_pool_get_array(ponds, "pondFreshWaterFlux", pondFreshWaterFlux)
       call MPAS_pool_get_array(ponds, "pondSnowDepthDifference", pondSnowDepthDifference)
       call MPAS_pool_get_array(ponds, "pondLidMeltFluxFraction", pondLidMeltFluxFraction)

       call MPAS_pool_get_array(diagnostics, "meltOnset", meltOnset)
       call MPAS_pool_get_array(diagnostics, "freezeOnset", freezeOnset)
       call MPAS_pool_get_array(diagnostics, "snowIceInterfaceTemperature", snowIceInterfaceTemperature)

       call MPAS_pool_get_array(snow, "snowLossToLeads", snowLossToLeads)
       call MPAS_pool_get_array(snow, "snowMeltMassCell", snowMeltMassCell)
       call MPAS_pool_get_array(snow, "snowMeltMassCategory", snowMeltMassCategory)

       ! high frequency coupling needs to cell center velocity
       if (config_use_high_frequency_coupling) then
          call seaice_interpolate_vertex_to_cell(mesh, boundary, uVelocityCell, uVelocity)
          call seaice_interpolate_vertex_to_cell(mesh, boundary, vVelocityCell, vVelocity)
       endif

       ! aerosols
       if (config_use_aerosols) then

          allocate(specificSnowAerosol(nAerosols, 2, nCategories))
          allocate(specificIceAerosol(nAerosols, 2, nCategories))

       else

          allocate(specificSnowAerosol(1, 1, 1))
          allocate(specificIceAerosol(1, 1, 1))
          specificSnowAerosol = 0.0_RKIND
          specificIceAerosol  = 0.0_RKIND

       endif

       ! code abort
       !abortFlag = .false.
       !abortMessage = ""

       !$offomp parallel do default(shared) private(iCategory,iAerosol,northernHemisphereMask,&
       !$offomp&  abortMessage) firstprivate(specificSnowAerosol,specificIceAerosol)
       !$offomp&  reduction(.or.:abortFlag)
       do iCell = 1, nCellsSolve

          ! code abort
          abortFlag = .false.
          abortMessage = ""

          ! initial state values
          iceAreaCellInitial(iCell) = iceAreaCell(iCell)

          do iCategory = 1, nCategories

             iceAreaCategoryInitial(iCategory,iCell) = iceAreaCategory(1,iCategory,iCell)
             iceVolumeCategoryInitial(iCategory,iCell) = iceVolumeCategory(1,iCategory,iCell)
             snowVolumeCategoryInitial(iCategory,iCell) = snowVolumeCategory(1,iCategory,iCell)

          enddo ! iCategory

          ! aerosol
          if (config_use_aerosols) then

             do iCategory = 1, nCategories
                do iAerosol = 1, nAerosols

                   specificSnowAerosol(iAerosol, 1, iCategory) = &
                        snowScatteringAerosol(iAerosol,iCategory,iCell) * snowVolumeCategoryInitial(iCategory,iCell)
                   specificSnowAerosol(iAerosol, 2, iCategory) = &
                        snowBodyAerosol(iAerosol,iCategory,iCell)       * snowVolumeCategoryInitial(iCategory,iCell)

                   specificIceAerosol(iAerosol, 1, iCategory) = &
                        iceScatteringAerosol(iAerosol,iCategory,iCell)   * iceVolumeCategoryInitial(iCategory,iCell)
                   specificIceAerosol(iAerosol, 2, iCategory) = &
                        iceBodyAerosol(iAerosol,iCategory,iCell)         * iceVolumeCategoryInitial(iCategory,iCell)

                enddo ! iAerosol
             enddo ! iCategory

          end if

          ! hemisphere mask
          if (latCell(iCell) > 0.0_RKIND) then
             northernHemisphereMask = .true.
          else
             northernHemisphereMask = .false.
          endif

          call icepack_warnings_clear()
          call icepack_step_therm1(&
               dt=config_dt, &
               aicen_init=iceAreaCategoryInitial(:,iCell), &
               vicen_init=iceVolumeCategoryInitial(:,iCell), &
               vsnon_init=snowVolumeCategoryInitial(:,iCell), &
               aice=iceAreaCell(iCell), &
               aicen=iceAreaCategory(1,:,iCell), &
               vice=iceVolumeCell(iCell), &
               vicen=iceVolumeCategory(1,:,iCell), &
               vsno=snowVolumeCell(iCell), &
               vsnon=snowVolumeCategory(1,:,iCell), &
               uvel=uVelocityCell(iCell), &
               vvel=vVelocityCell(iCell), &
               Tsfc=surfaceTemperature(1,:,iCell), &
               zqsn=snowEnthalpy(:,:,iCell), &
               zqin=iceEnthalpy(:,:,iCell), &
               zSin=iceSalinity(:,:,iCell), &
               alvl=levelIceArea(1,:,iCell), &
               vlvl=levelIceVolume(1,:,iCell), &
               apnd=pondArea(1,:,iCell), &
               hpnd=pondDepth(1,:,iCell), &
               ipnd=pondLidThickness(1,:,iCell), &
               iage=iceAge(1,:,iCell), &
               FY=firstYearIceArea(1,:,iCell), &
               aerosno=specificSnowAerosol(:,:,:), &
               aeroice=specificIceAerosol(:,:,:), &
               uatm=uAirVelocity(iCell), &
               vatm=vAirVelocity(iCell), &
               wind=windSpeed(iCell), &
               zlvl=airLevelHeight(iCell), &
               Qa=airSpecificHumidity(iCell), &
               rhoa=airDensity(iCell), &
               Tair=airTemperature(iCell), &
               Tref=atmosReferenceTemperature2m(iCell), &
               Qref=atmosReferenceHumidity2m(iCell), &
               Uref=atmosReferenceSpeed10m(iCell), &
               Cdn_atm_ratio=airOceanDragCoefficientRatio(iCell), &
               Cdn_ocn=oceanDragCoefficient(iCell), &
               Cdn_ocn_skin=oceanDragCoefficientSkin(iCell), &
               Cdn_ocn_floe=oceanDragCoefficientFloe(iCell), &
               Cdn_ocn_keel=oceanDragCoefficientKeel(iCell), &
               Cdn_atm=airDragCoefficient(iCell), &
               Cdn_atm_skin=airDragCoefficientSkin(iCell), &
               Cdn_atm_floe=airDragCoefficientFloe(iCell), &
               Cdn_atm_pond=airDragCoefficientPond(iCell), &
               Cdn_atm_rdg=airDragCoefficientRidge(iCell), &
               hfreebd=dragFreeboard(iCell), &
               hdraft=dragIceSnowDraft(iCell), &
               hridge=dragRidgeHeight(iCell), &
               distrdg=dragRidgeSeparation(iCell), &
               hkeel=dragKeelDepth(iCell), &
               dkeel=dragKeelSeparation(iCell), &
               lfloe=dragFloeLength(iCell), &
               dfloe=dragFloeSeparation(iCell), &
               strax=airStressForcingU(iCell), &
               stray=airStressForcingV(iCell), &
               strairxT=airStressCellU(iCell), &
               strairyT=airStressCellV(iCell), &
               potT=airPotentialTemperature(iCell), &
               sst=seaSurfaceTemperature(iCell), &
               sss=seaSurfaceSalinity(iCell), &
               Tf=seaFreezingTemperature(iCell), &
               strocnxT=oceanStressCellU(iCell), &
               strocnyT=oceanStressCellV(iCell), &
               fbot=oceanHeatFluxIceBottom(iCell), &
               Tbot=bottomTemperature(iCell), &
               Tsnice=snowIceInterfaceTemperature(iCell), &
               frzmlt=freezingMeltingPotential(iCell), &
               rsiden=lateralIceMeltFractionCategory(:,iCell), &
               wlat=lateralIceMeltRate(iCell), &
               fsnow=snowfallRate(iCell), &
               frain=rainfallRate(iCell), &
               fpond=pondFreshWaterFlux(iCell), &
               fsloss=snowLossToLeads(iCell), &
               fsurf=surfaceHeatFlux(iCell), &
               fsurfn=surfaceHeatFluxCategory(:,iCell), &
               fcondtop=surfaceConductiveFlux(iCell), &
               fcondtopn=surfaceConductiveFluxCategory(:,iCell), &
               fcondbot=bottomConductiveFlux(iCell), &
               fcondbotn=bottomConductiveFluxCategory(:,iCell), &
               fswsfcn=surfaceShortwaveFlux(:,iCell), &
               fswintn=interiorShortwaveFlux(:,iCell), &
               fswthrun=penetratingShortwaveFlux(:,iCell), &
               fswabs=absorbedShortwaveFlux(iCell), &
               flwout=longwaveUp(iCell), &
               Sswabsn=absorbedShortwaveSnowLayer(:,:,iCell), &
               Iswabsn=absorbedShortwaveIceLayer(:,:,iCell), &
               flw=longwaveDown(iCell), &
               fsens=sensibleHeatFlux(iCell), &
               fsensn=sensibleHeatFluxCategory(:,iCell), &
               flat=latentHeatFlux(iCell), &
               flatn=latentHeatFluxCategory(:,iCell), &
               evap=evaporativeWaterFlux(iCell), &
               evaps=evaporativeWaterFluxSnow(iCell), &
               evapi=evaporativeWaterFluxIce(iCell), &
               fresh=oceanFreshWaterFlux(iCell), &
               fsalt=oceanSaltFlux(iCell), &
               fhocn=oceanHeatFlux(iCell), &
               fswthru=oceanShortwaveFlux(iCell), & !??
               flatn_f=latentHeatFluxCouple(:,iCell), &
               fsensn_f=sensibleHeatFluxCouple(:,iCell), &
               fsurfn_f=surfaceHeatFluxCouple(:,iCell), &
               fcondtopn_f=surfaceConductiveFluxCouple(:,iCell) , &
               faero_atm=atmosAerosolFlux(:,iCell), &
               faero_ocn=oceanAerosolFlux(:,iCell), &
               dhsn=pondSnowDepthDifference(:,iCell), &
               ffracn=pondLidMeltFluxFraction(:,iCell), &
               meltt=surfaceIceMelt(iCell), &
               melttn=surfaceIceMeltCategory(:,iCell), &
               meltb=basalIceMelt(iCell), &
               meltbn=basalIceMeltCategory(:,iCell), &
               melts=snowMelt(iCell), &
               meltsn=snowMeltCategory(:,iCell), &
               congel=congelation(iCell), &
               congeln=congelationCategory(:,iCell), &
               snoice=snowiceFormation(iCell), &
               snoicen=snowiceFormationCategory(:,iCell), &
               dsnow=snowThicknessChangeCell(iCell), &
               dsnown=snowThicknessChangeCategory(:,iCell), &
               meltsliq=snowMeltMassCell(iCell), &
               meltsliqn=snowMeltMassCategory(:,iCell), &
               rsnwn=snowGrainRadius(:,:,iCell), &
               smicen=snowIceMass(:,:,iCell), &
               smliqn=snowLiquidMass(:,:,iCell), &
               lmask_n=northernHemisphereMask, &
               lmask_s=(.not. northernHemisphereMask), &
               mlt_onset=meltOnset(iCell), &
               frz_onset=freezeOnset(iCell), &
               yday=dayOfYear, &
               prescribed_ice=config_use_prescribed_ice)

          abortFlag = icepack_warnings_aborted()
          call seaice_icepack_write_warnings(abortFlag)

          ! cell-specific abort message
          if (abortFlag) then
             call mpas_log_write("column_vertical_thermodynamics: "//trim(abortMessage) , messageType=MPAS_LOG_ERR)
             call mpas_log_write("iCell: $i", messageType=MPAS_LOG_ERR, intArgs=(/indexToCellID(iCell)/))

             call mpas_log_write("config_dt: $r", messageType=MPAS_LOG_ERR, realArgs=(/config_dt/))
             call mpas_log_write("nCategories: $i", messageType=MPAS_LOG_ERR, intArgs=(/nCategories/))
             call mpas_log_write("nIceLayers: $i", messageType=MPAS_LOG_ERR, intArgs=(/nIceLayers/))
             call mpas_log_write("nSnowLayers: $i", messageType=MPAS_LOG_ERR, intArgs=(/nSnowLayers/))
             call mpas_log_write("nAerosols: $i", messageType=MPAS_LOG_ERR, intArgs=(/nAerosols/))
             call mpas_log_write("openWaterArea: $r", messageType=MPAS_LOG_ERR, realArgs=(/openWaterArea(iCell)/))
             call mpas_log_write("iceAreaCategoryInitial: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/iceAreaCategoryInitial(:,iCell)/))
             call mpas_log_write("iceVolumeCategoryInitial: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/iceVolumeCategoryInitial(:,iCell)/))
             call mpas_log_write("snowVolumeCategoryInitial: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/snowVolumeCategoryInitial(:,iCell)/))
             call mpas_log_write("iceAreaCell: $r", messageType=MPAS_LOG_ERR, realArgs=(/iceAreaCell(iCell)/))
             call mpas_log_write("iceAreaCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/iceAreaCategory(1,:,iCell)/))
             call mpas_log_write("iceVolumeCell: $r", messageType=MPAS_LOG_ERR, realArgs=(/iceVolumeCell(iCell)/))
             call mpas_log_write("iceVolumeCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/iceVolumeCategory(1,:,iCell)/))
             call mpas_log_write("snowVolumeCell: $r", messageType=MPAS_LOG_ERR, realArgs=(/snowVolumeCell(iCell)/))
             call mpas_log_write("snowVolumeCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/snowVolumeCategory(1,:,iCell)/))
             call mpas_log_write("uVelocityCell: $r", messageType=MPAS_LOG_ERR, realArgs=(/uVelocityCell(iCell)/))
             call mpas_log_write("vVelocityCell: $r", messageType=MPAS_LOG_ERR, realArgs=(/vVelocityCell(iCell)/))
             call mpas_log_write("surfaceTemperature: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/surfaceTemperature(1,:,iCell)/))
             do iCategory = 1, nCategories
                call mpas_log_write("snowEnthalpy: $i "//repeat("$r ", nSnowLayers), messageType=MPAS_LOG_ERR, intArgs=(/iCategory/), realArgs=(/snowEnthalpy(:,iCategory,iCell)/))
             enddo ! iCategory
             do iCategory = 1, nCategories
                call mpas_log_write("iceEnthalpy: $i "//repeat("$r ", nIceLayers), messageType=MPAS_LOG_ERR, intArgs=(/iCategory/), realArgs=(/iceEnthalpy(:,iCategory,iCell)/))
             enddo ! iCategory
             do iCategory = 1, nCategories
                call mpas_log_write("iceSalinity: $i "//repeat("$r ", nIceLayers), messageType=MPAS_LOG_ERR, intArgs=(/iCategory/), realArgs=(/iceSalinity(:,iCategory,iCell)/))
             enddo ! iCategory
             do iCategory = 1, nCategories
                call mpas_log_write("snowIceMass: $i "//repeat("$r ", nSnowLayers), messageType=MPAS_LOG_ERR, intArgs=(/iCategory/), realArgs=(/snowIceMass(:,iCategory,iCell)/))
             enddo ! iCategory
             do iCategory = 1, nCategories
                call mpas_log_write("snowLiquidMass: $i "//repeat("$r ", nSnowLayers), messageType=MPAS_LOG_ERR, intArgs=(/iCategory/), realArgs=(/snowLiquidMass(:,iCategory,iCell)/))
             enddo ! iCategory
             call mpas_log_write("levelIceArea: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/levelIceArea(1,:,iCell)/))
             call mpas_log_write("levelIceVolume: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/levelIceVolume(1,:,iCell)/))
             call mpas_log_write("pondArea: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/pondArea(1,:,iCell)/))
             call mpas_log_write("pondDepth: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/pondDepth(1,:,iCell)/))
             call mpas_log_write("pondLidThickness: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/pondLidThickness(1,:,iCell)/))
             call mpas_log_write("iceAge: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/iceAge(1,:,iCell)/))
             call mpas_log_write("firstYearIceArea: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/firstYearIceArea(1,:,iCell)/))
             do iCategory = 1, nCategories
                call mpas_log_write("snowGrainRadius: $i "//repeat("$r ", nSnowLayers), messageType=MPAS_LOG_ERR, intArgs=(/iCategory/), realArgs=(/snowGrainRadius(:,iCategory,iCell)/))
             enddo ! iCategory
             call mpas_log_write("config_use_snow_liquid_ponds: $l", messageType=MPAS_LOG_ERR, logicArgs=(/config_use_snow_liquid_ponds/))
             if (config_use_aerosols) then
                do iCategory = 1, nCategories
                   call mpas_log_write("specificSnowAerosol $i 1: "//repeat("$r ", nAerosols), messageType=MPAS_LOG_ERR, intArgs=(/iCategory/), realArgs=(/specificSnowAerosol(:,1,iCategory)/))
                   call mpas_log_write("specificSnowAerosol $i 2: "//repeat("$r ", nAerosols), messageType=MPAS_LOG_ERR, intArgs=(/iCategory/), realArgs=(/specificSnowAerosol(:,2,iCategory)/))
                enddo ! iCategory
                do iCategory = 1, nCategories
                   call mpas_log_write("specificIceAerosol $i 1: "//repeat("$r ", nAerosols), messageType=MPAS_LOG_ERR, intArgs=(/iCategory/), realArgs=(/specificIceAerosol(:,1,iCategory)/))
                   call mpas_log_write("specificIceAerosol $i 2: "//repeat("$r ", nAerosols), messageType=MPAS_LOG_ERR, intArgs=(/iCategory/), realArgs=(/specificIceAerosol(:,2,iCategory)/))
                enddo ! iCategory
             endif
             call mpas_log_write("uAirVelocity: $r", messageType=MPAS_LOG_ERR, realArgs=(/uAirVelocity(iCell)/))
             call mpas_log_write("vAirVelocity: $r", messageType=MPAS_LOG_ERR, realArgs=(/vAirVelocity(iCell)/))
             call mpas_log_write("windSpeed: $r", messageType=MPAS_LOG_ERR, realArgs=(/windSpeed(iCell)/))
             call mpas_log_write("airLevelHeight: $r", messageType=MPAS_LOG_ERR, realArgs=(/airLevelHeight(iCell)/))
             call mpas_log_write("airSpecificHumidity: $r", messageType=MPAS_LOG_ERR, realArgs=(/airSpecificHumidity(iCell)/))
             call mpas_log_write("airDensity: $r", messageType=MPAS_LOG_ERR, realArgs=(/airDensity(iCell)/))
             call mpas_log_write("airTemperature: $r", messageType=MPAS_LOG_ERR, realArgs=(/airTemperature(iCell)/))
             call mpas_log_write("atmosReferenceTemperature2m: $r", messageType=MPAS_LOG_ERR, realArgs=(/atmosReferenceTemperature2m(iCell)/))
             call mpas_log_write("atmosReferenceHumidity2m: $r", messageType=MPAS_LOG_ERR, realArgs=(/atmosReferenceHumidity2m(iCell)/))
             call mpas_log_write("atmosReferenceSpeed10m: $r", messageType=MPAS_LOG_ERR, realArgs=(/atmosReferenceSpeed10m(iCell)/))
             call mpas_log_write("airOceanDragCoefficientRatio: $r", messageType=MPAS_LOG_ERR, realArgs=(/airOceanDragCoefficientRatio(iCell)/))
             call mpas_log_write("oceanDragCoefficient: $r", messageType=MPAS_LOG_ERR, realArgs=(/oceanDragCoefficient(iCell)/))
             call mpas_log_write("oceanDragCoefficientSkin: $r", messageType=MPAS_LOG_ERR, realArgs=(/oceanDragCoefficientSkin(iCell)/))
             call mpas_log_write("oceanDragCoefficientFloe: $r", messageType=MPAS_LOG_ERR, realArgs=(/oceanDragCoefficientFloe(iCell)/))
             call mpas_log_write("oceanDragCoefficientKeel: $r", messageType=MPAS_LOG_ERR, realArgs=(/oceanDragCoefficientKeel(iCell)/))
             call mpas_log_write("airDragCoefficient: $r", messageType=MPAS_LOG_ERR, realArgs=(/airDragCoefficient(iCell)/))
             call mpas_log_write("airDragCoefficientSkin: $r", messageType=MPAS_LOG_ERR, realArgs=(/airDragCoefficientSkin(iCell)/))
             call mpas_log_write("airDragCoefficientFloe: $r", messageType=MPAS_LOG_ERR, realArgs=(/airDragCoefficientFloe(iCell)/))
             call mpas_log_write("airDragCoefficientPond: $r", messageType=MPAS_LOG_ERR, realArgs=(/airDragCoefficientPond(iCell)/))
             call mpas_log_write("airDragCoefficientRidge: $r", messageType=MPAS_LOG_ERR, realArgs=(/airDragCoefficientRidge(iCell)/))
             call mpas_log_write("dragFreeboard: $r", messageType=MPAS_LOG_ERR, realArgs=(/dragFreeboard(iCell)/))
             call mpas_log_write("dragIceSnowDraft: $r", messageType=MPAS_LOG_ERR, realArgs=(/dragIceSnowDraft(iCell)/))
             call mpas_log_write("dragRidgeHeight: $r", messageType=MPAS_LOG_ERR, realArgs=(/dragRidgeHeight(iCell)/))
             call mpas_log_write("dragRidgeSeparation: $r", messageType=MPAS_LOG_ERR, realArgs=(/dragRidgeSeparation(iCell)/))
             call mpas_log_write("dragKeelDepth: $r", messageType=MPAS_LOG_ERR, realArgs=(/dragKeelDepth(iCell)/))
             call mpas_log_write("dragKeelSeparation: $r", messageType=MPAS_LOG_ERR, realArgs=(/dragKeelSeparation(iCell)/))
             call mpas_log_write("dragFloeLength: $r", messageType=MPAS_LOG_ERR, realArgs=(/dragFloeLength(iCell)/))
             call mpas_log_write("dragFloeSeparation: $r", messageType=MPAS_LOG_ERR, realArgs=(/dragFloeSeparation(iCell)/))
             call mpas_log_write("airStressForcingU: $r", messageType=MPAS_LOG_ERR, realArgs=(/airStressForcingU(iCell)/))
             call mpas_log_write("airStressForcingV: $r", messageType=MPAS_LOG_ERR, realArgs=(/airStressForcingV(iCell)/))
             call mpas_log_write("airStressCellU: $r", messageType=MPAS_LOG_ERR, realArgs=(/airStressCellU(iCell)/))
             call mpas_log_write("airStressCellV: $r", messageType=MPAS_LOG_ERR, realArgs=(/airStressCellV(iCell)/))
             call mpas_log_write("airPotentialTemperature: $r", messageType=MPAS_LOG_ERR, realArgs=(/airPotentialTemperature(iCell)/))
             call mpas_log_write("seaSurfaceTemperature: $r", messageType=MPAS_LOG_ERR, realArgs=(/seaSurfaceTemperature(iCell)/))
             call mpas_log_write("seaSurfaceSalinity: $r", messageType=MPAS_LOG_ERR, realArgs=(/seaSurfaceSalinity(iCell)/))
             call mpas_log_write("seaFreezingTemperature: $r", messageType=MPAS_LOG_ERR, realArgs=(/seaFreezingTemperature(iCell)/))
             call mpas_log_write("oceanStressCellU: $r", messageType=MPAS_LOG_ERR, realArgs=(/oceanStressCellU(iCell)/))
             call mpas_log_write("oceanStressCellV: $r", messageType=MPAS_LOG_ERR, realArgs=(/oceanStressCellV(iCell)/))
             call mpas_log_write("oceanHeatFluxIceBottom: $r", messageType=MPAS_LOG_ERR, realArgs=(/oceanHeatFluxIceBottom(iCell)/))
             call mpas_log_write("freezingMeltingPotential: $r", messageType=MPAS_LOG_ERR, realArgs=(/freezingMeltingPotential(iCell)/))
             call mpas_log_write("lateralIceMeltFractionCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/lateralIceMeltFractionCategory(:,iCell)/))
             call mpas_log_write("lateralIceMeltRate: $r", messageType=MPAS_LOG_ERR, realArgs=(/lateralIceMeltRate(iCell)/))
             call mpas_log_write("snowfallRate: $r", messageType=MPAS_LOG_ERR, realArgs=(/snowfallRate(iCell)/))
             call mpas_log_write("rainfallRate: $r", messageType=MPAS_LOG_ERR, realArgs=(/rainfallRate(iCell)/))
             call mpas_log_write("pondFreshWaterFlux: $r", messageType=MPAS_LOG_ERR, realArgs=(/pondFreshWaterFlux(iCell)/))
             call mpas_log_write("snowLossToLeads: $r", messageType=MPAS_LOG_ERR, realArgs=(/snowLossToLeads(iCell)/))
             call mpas_log_write("surfaceHeatFlux: $r", messageType=MPAS_LOG_ERR, realArgs=(/surfaceHeatFlux(iCell)/))
             call mpas_log_write("surfaceHeatFluxCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/surfaceHeatFluxCategory(:,iCell)/))
             call mpas_log_write("surfaceConductiveFlux: $r", messageType=MPAS_LOG_ERR, realArgs=(/surfaceConductiveFlux(iCell)/))
             call mpas_log_write("surfaceConductiveFluxCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/surfaceConductiveFluxCategory(:,iCell)/))
             call mpas_log_write("surfaceShortwaveFlux: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/surfaceShortwaveFlux(:,iCell)/))
             call mpas_log_write("interiorShortwaveFlux: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/interiorShortwaveFlux(:,iCell)/))
             call mpas_log_write("penetratingShortwaveFlux: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/penetratingShortwaveFlux(:,iCell)/))
             call mpas_log_write("absorbedShortwaveFlux: $r", messageType=MPAS_LOG_ERR, realArgs=(/absorbedShortwaveFlux(iCell)/))
             call mpas_log_write("longwaveUp: $r", messageType=MPAS_LOG_ERR, realArgs=(/longwaveUp(iCell)/))
             do iCategory = 1, nCategories
                call mpas_log_write("absorbedShortwaveSnowLayer: $i "//repeat("$r ", nSnowLayers), messageType=MPAS_LOG_ERR, intArgs=(/iCategory/), realArgs=(/absorbedShortwaveSnowLayer(:,iCategory,iCell)/))
             enddo ! iCategory
             do iCategory = 1, nCategories
                call mpas_log_write("absorbedShortwaveIceLayer: $i "//repeat("$r ", nIceLayers), messageType=MPAS_LOG_ERR, intArgs=(/iCategory/), realArgs=(/absorbedShortwaveIceLayer(:,iCategory,iCell)/))
             enddo ! iCategory
             call mpas_log_write("longwaveDown: $r", messageType=MPAS_LOG_ERR, realArgs=(/longwaveDown(iCell)/))
             call mpas_log_write("solarZenithAngleCosine: $r", messageType=MPAS_LOG_ERR, realArgs=(/solarZenithAngleCosine(iCell)/))
             call mpas_log_write("sensibleHeatFlux: $r", messageType=MPAS_LOG_ERR, realArgs=(/sensibleHeatFlux(iCell)/))
             call mpas_log_write("sensibleHeatFluxCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/sensibleHeatFluxCategory(:,iCell)/))
             call mpas_log_write("latentHeatFlux: $r", messageType=MPAS_LOG_ERR, realArgs=(/latentHeatFlux(iCell)/))
             call mpas_log_write("latentHeatFluxCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/latentHeatFluxCategory(:,iCell)/))
             call mpas_log_write("evaporativeWaterFlux: $r", messageType=MPAS_LOG_ERR, realArgs=(/evaporativeWaterFlux(iCell)/))
             call mpas_log_write("oceanFreshWaterFlux: $r", messageType=MPAS_LOG_ERR, realArgs=(/oceanFreshWaterFlux(iCell)/))
             call mpas_log_write("oceanSaltFlux: $r", messageType=MPAS_LOG_ERR, realArgs=(/oceanSaltFlux(iCell)/))
             call mpas_log_write("oceanHeatFlux: $r", messageType=MPAS_LOG_ERR, realArgs=(/oceanHeatFlux(iCell)/))
             call mpas_log_write("oceanShortwaveFlux: $r", messageType=MPAS_LOG_ERR, realArgs=(/oceanShortwaveFlux(iCell)/))
             call mpas_log_write("latentHeatFluxCouple: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/latentHeatFluxCouple(:,iCell)/))
             call mpas_log_write("sensibleHeatFluxCouple: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/sensibleHeatFluxCouple(:,iCell)/))
             call mpas_log_write("surfaceHeatFluxCouple: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/surfaceHeatFluxCouple(:,iCell)/))
             call mpas_log_write("surfaceConductiveFluxCouple: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/surfaceConductiveFluxCouple(:,iCell)/))
             call mpas_log_write("atmosAerosolFlux: "//repeat("$r ", nAerosols), messageType=MPAS_LOG_ERR, realArgs=(/atmosAerosolFlux(:,iCell)/))
             call mpas_log_write("oceanAerosolFlux: "//repeat("$r ", nAerosols), messageType=MPAS_LOG_ERR, realArgs=(/oceanAerosolFlux(:,iCell)/))
             call mpas_log_write("pondSnowDepthDifference: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/pondSnowDepthDifference(:,iCell)/))
             call mpas_log_write("pondLidMeltFluxFraction: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/pondLidMeltFluxFraction(:,iCell)/))
             call mpas_log_write("surfaceIceMelt: $r", messageType=MPAS_LOG_ERR, realArgs=(/surfaceIceMelt(iCell)/))
             call mpas_log_write("surfaceIceMeltCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/surfaceIceMeltCategory(:,iCell)/))
             call mpas_log_write("basalIceMelt: $r", messageType=MPAS_LOG_ERR, realArgs=(/basalIceMelt(iCell)/))
             call mpas_log_write("basalIceMeltCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/basalIceMeltCategory(:,iCell)/))
             call mpas_log_write("lateralIceMelt: $r", messageType=MPAS_LOG_ERR, realArgs=(/lateralIceMelt(iCell)/))
             call mpas_log_write("snowMelt: $r", messageType=MPAS_LOG_ERR, realArgs=(/snowMelt(iCell)/))
             call mpas_log_write("snowMeltCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/snowMeltCategory(:,iCell)/))
             call mpas_log_write("snowMeltMassCell: $r", messageType=MPAS_LOG_ERR, realArgs=(/snowMeltMassCell(iCell)/))
             call mpas_log_write("snowMeltMassCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/snowMeltMassCategory(:,iCell)/))
             call mpas_log_write("congelation: $r", messageType=MPAS_LOG_ERR, realArgs=(/congelation(iCell)/))
             call mpas_log_write("congelationCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/congelationCategory(:,iCell)/))
             call mpas_log_write("snowiceFormation: $r", messageType=MPAS_LOG_ERR, realArgs=(/snowiceFormation(iCell)/))
             call mpas_log_write("snowiceFormationCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/snowiceFormationCategory(:,iCell)/))
             call mpas_log_write("snowThicknessChangeCategory: "//repeat("$r ", nCategories), messageType=MPAS_LOG_ERR, realArgs=(/snowThicknessChangeCategory(:,iCell)/))
             call mpas_log_write("frazilFormation: $r", messageType=MPAS_LOG_ERR, realArgs=(/frazilFormation(iCell)/))
             call mpas_log_write("northernHemisphereMask: $l", messageType=MPAS_LOG_ERR, logicArgs=(/northernHemisphereMask/))
             call mpas_log_write("meltOnset: $r", messageType=MPAS_LOG_ERR, realArgs=(/meltOnset(iCell)/))
             call mpas_log_write("freezeOnset: $r", messageType=MPAS_LOG_ERR, realArgs=(/freezeOnset(iCell)/))
             call mpas_log_write("dayOfYear: $r", messageType=MPAS_LOG_ERR, realArgs=(/dayOfYear/))
             call mpas_log_write("config_use_prescribed_ice: $l", messageType=MPAS_LOG_ERR, logicArgs=(/config_use_prescribed_ice/))
          endif

          ! aerosol
          if (config_use_aerosols) then

             do iCategory = 1, nCategories
                do iAerosol = 1, nAerosols

                   if (snowVolumeCategory(1,iCategory,iCell) > seaicePuny) &
                        specificSnowAerosol(iAerosol, :, iCategory) = &
                        specificSnowAerosol(iAerosol, :, iCategory) / snowVolumeCategory(1,iCategory,iCell)

                   if (iceVolumeCategory(1,iCategory,iCell) > seaicePuny) &
                        specificIceAerosol(iAerosol, :, iCategory)  = &
                        specificIceAerosol(iAerosol, :, iCategory)  / iceVolumeCategory(1,iCategory,iCell)

                   snowScatteringAerosol(iAerosol,iCategory,iCell) = specificSnowAerosol(iAerosol, 1, iCategory)
                   snowBodyAerosol(iAerosol,iCategory,iCell)       = specificSnowAerosol(iAerosol, 2, iCategory)

                   iceScatteringAerosol(iAerosol,iCategory,iCell)  = specificIceAerosol(iAerosol, 1, iCategory)
                   iceBodyAerosol(iAerosol,iCategory,iCell)        = specificIceAerosol(iAerosol, 2, iCategory)

                enddo ! iAerosol
             enddo ! iCategory

          endif
       enddo ! iCell

       ! error-checking
       call seaice_critical_error_write_block(domain, block, abortFlag)
       call seaice_check_critical_error(domain, abortFlag)

       ! aerosols
       deallocate(specificSnowAerosol)
       deallocate(specificIceAerosol)

       block => block % next
    end do

  end subroutine column_vertical_thermodynamics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_itd_thermodynamics
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine column_itd_thermodynamics(domain, clock)

    use icepack_intfc, only: &
         icepack_step_therm2

    use seaice_constants, only: &
         seaicePuny

    type(domain_type), intent(inout) :: domain

    type(MPAS_clock_type), intent(in) :: clock

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         icestate, &
         tracers, &
         tracers_aggregate, &
         atmos_coupling, &
         ocean_coupling, &
!         wave_coupling, &
!         floe_size_distribution, &
         ocean_fluxes, &
         melt_growth_rates, &
         ponds, &
         biogeochemistry, &
         initial, &
         diagnostics, &
         aerosols

    ! configs
    real(kind=RKIND), pointer :: &
         config_dt

    logical, pointer :: &
         config_update_ocean_fluxes, &
         config_use_column_biogeochemistry, &
         config_use_zaerosols

    ! dimensions
    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
!         nFloeCategories, &
         nIceLayers, &
         nSnowLayers, &
         nAerosols, &
         nBioLayers, &
         nBioLayersP1

    ! variables
    real(kind=RKIND), dimension(:), pointer :: &
         openWaterArea, &
         iceAreaCell, &
         seaFreezingTemperature, &
         seaSurfaceSalinity, &
         lateralIceMeltRate, &
         lateralIceMelt, &
         freezingMeltingPotential, &
         frazilFormation, &
         rainfallRate, &
         pondFreshWaterFlux, &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         freezeOnset, &
         categoryThicknessLimits, &
         frazilGrowthDiagnostic
!         frazilGrowthDiagnostic, &
!         WaveFrequency,      &
!         WaveFreqBinWidth,   &
!         SignificantWaveHeight,&
!         FloeSizeBinCenter,  &
!         FloeSizeBinWidth

    real(kind=RKIND), dimension(:,:), pointer :: &
         iceAreaCategoryInitial, &
         iceVolumeCategoryInitial, &
         lateralIceMeltFractionCategory, &
         oceanAerosolFlux, &
         oceanBioFluxes, &
         oceanBioConcentrations, &
         oceanBioConcentrationsInUse, &
         initialSalinityProfile
!         WaveSpectra, &
!         DFloeSizeNewIce, &
!         DFloeSizeLateralGrowth, &
!         DFloeSizeLateralMelt, &
!         DFloeSizeWaves, &
!         DFloeSizeWeld

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         brineFraction

    integer, dimension(:,:), pointer :: &
         newlyFormedIce

    integer, dimension(:), pointer :: &
         indexToCellID

    ! local
    integer :: &
         iCell, &
         iCategory, &
         iBioTracers, &
         iBioData, &
         iBioLayers

    real(kind=RKIND), dimension(:,:), allocatable :: &
         oceanBioFluxesTemp

    logical, dimension(:), allocatable :: &
         newlyFormedIceLogical

    logical :: &
         abortFlag, &
         setGetPhysicsTracers, &
         setGetBGCTracers

    character(len=strKIND) :: &
         abortMessage, &
         abortLocation

    real(kind=RKIND) :: &
         dayOfYear

    ! day of year
    call get_day_of_year(clock, dayOfYear)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", ocean_fluxes)
       call MPAS_pool_get_subpool(block % structs, "melt_growth_rates", melt_growth_rates)
       call MPAS_pool_get_subpool(block % structs, "ponds", ponds)
       call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)
       call MPAS_pool_get_subpool(block % structs, "initial", initial)
       call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnostics)
       call MPAS_pool_get_subpool(block % structs, "aerosols", aerosols)
!       call MPAS_pool_get_subpool(block % structs, "wave_coupling", wave_coupling)
!       call MPAS_pool_get_subpool(block % structs, "floe_size_distribution", floe_size_distribution)

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_update_ocean_fluxes", config_update_ocean_fluxes)
       call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
       call MPAS_pool_get_config(block % configs, "config_use_zaerosols", config_use_zaerosols)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
!       call MPAS_pool_get_dimension(mesh, "nFloeCategories", nFloeCategories)
       call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)
       call MPAS_pool_get_dimension(mesh, "nAerosols", nAerosols)
       call MPAS_pool_get_dimension(block % dimensions, "nBioLayers", nBioLayers)

       call MPAS_pool_get_dimension(block % dimensions, "nBioLayersP1", nBioLayersP1)
       call MPAS_pool_get_array(mesh, "indexToCellID", indexToCellID)

       call MPAS_pool_get_array(icestate, "iceAreaCategoryInitial", iceAreaCategoryInitial)
       call MPAS_pool_get_array(icestate, "iceVolumeCategoryInitial", iceVolumeCategoryInitial)
       call MPAS_pool_get_array(icestate, "openWaterArea", openWaterArea)

       call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "brineFraction", brineFraction, 1)

       call MPAS_pool_get_array(atmos_coupling, "rainfallRate", rainfallRate)

       call MPAS_pool_get_array(ocean_coupling, "freezingMeltingPotential", freezingMeltingPotential)
       call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)
       call MPAS_pool_get_array(ocean_coupling, "seaSurfaceSalinity", seaSurfaceSalinity)

!       call MPAS_pool_get_array(wave_coupling, "SignificantWaveHeight", SignificantWaveHeight)
!       call MPAS_pool_get_array(wave_coupling, "WaveSpectra", WaveSpectra)
!       call MPAS_pool_get_array(wave_coupling, "WaveFrequency", WaveFrequency)
!       call MPAS_pool_get_array(wave_coupling, "WaveFreqBinWidth", WaveFreqBinWidth)

!       call MPAS_pool_get_array(floe_size_distribution, "FloeSizeBinCenter", FloeSizeBinCenter)
!       call MPAS_pool_get_array(floe_size_distribution, "FloeSizeBinWidth", FloeSizeBinWidth)
!       call MPAS_pool_get_array(floe_size_distribution, "DFloeSizeNewIce", DFloeSizeNewIce)
!       call MPAS_pool_get_array(floe_size_distribution, "DFloeSizeLateralGrowth", DFloeSizeLateralGrowth)
!       call MPAS_pool_get_array(floe_size_distribution, "DFloeSizeLateralMelt", DFloeSizeLateralMelt)
!       call MPAS_pool_get_array(floe_size_distribution, "DFloeSizeWeld", DFloeSizeWeld)
!       call MPAS_pool_get_array(floe_size_distribution, "DFloeSizeWaves", DFloeSizeWaves)

       call MPAS_pool_get_array(ocean_fluxes, "oceanFreshWaterFlux", oceanFreshWaterFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanSaltFlux", oceanSaltFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanHeatFlux", oceanHeatFlux)

       call MPAS_pool_get_array(melt_growth_rates, "lateralIceMeltFractionCategory", lateralIceMeltFractionCategory)
       call MPAS_pool_get_array(melt_growth_rates, "lateralIceMeltRate", lateralIceMeltRate)
       call MPAS_pool_get_array(melt_growth_rates, "lateralIceMelt", lateralIceMelt)
       call MPAS_pool_get_array(melt_growth_rates, "frazilFormation", frazilFormation)
       call MPAS_pool_get_array(melt_growth_rates, "frazilGrowthDiagnostic", frazilGrowthDiagnostic)

       call MPAS_pool_get_array(ponds, "pondFreshWaterFlux", pondFreshWaterFlux)

       call MPAS_pool_get_array(aerosols, "oceanAerosolFlux", oceanAerosolFlux)

       call MPAS_pool_get_array(biogeochemistry, "newlyFormedIce", newlyFormedIce)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioFluxes", oceanBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioConcentrations", oceanBioConcentrations)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioConcentrationsInUse", oceanBioConcentrationsInUse)

       call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)
       call MPAS_pool_get_array(initial, "categoryThicknessLimits", categoryThicknessLimits)

       call MPAS_pool_get_array(diagnostics, "freezeOnset", freezeOnset)

       ! newly formed ice
       allocate(newlyFormedIceLogical(nCategories))
       allocate(oceanBioFluxesTemp(ciceTracerObject % nBioTracers,nCellsSolve))
       setGetPhysicsTracers = .true.
       setGetBGCTracers     = (config_use_column_biogeochemistry .or. config_use_zaerosols)

       oceanBioConcentrationsInUse(:,:) = 0.0_RKIND
       ! code abort
       abortFlag = .false.
       abortMessage = ""

       !$offomp parallel do default(shared) private(iCategory,iBioTracers,iBioData,&
       !$offomp&   totalCarbonInitial,abortMessage,oceanBioFluxesTemp,totalCarbonFinal,&
       !$offomp&   totalCarbonCatInitial, totalCarbonCatFinal, &
       !$offomp&   oceanCarbonFlux, carbonError) &
       !$offomp&   firstprivate(newlyFormedIceLogical) &
       !$offomp&   reduction(.or.:abortFlag)
       do iCell = 1, nCellsSolve

          ! newly formed ice
          do iCategory = 1, nCategories
             newlyFormedIceLogical(iCategory) = (newlyFormedIce(iCategory,iCell) == 1)
          enddo ! iCategory

          ! read the required ocean concentration fields into the allocated array
          do iBioTracers = 1, ciceTracerObject % nBioTracers
             iBioData = ciceTracerObject % index_LayerIndexToDataArray(iBioTracers)
             oceanBioConcentrationsInUse(iBioTracers, iCell) = oceanBioConcentrations(iBioData, iCell)
          enddo ! iBioTracers

          ! set the category tracer array
          call set_cice_tracer_array_category(block, ciceTracerObject,&
                 tracerArrayCategory, iCell, setGetPhysicsTracers, setGetBGCTracers)

          oceanBioFluxesTemp(:,iCell) = 0.0_RKIND

          ! CICE calculates the Sig Wave Height from the wave spectrum as follows before step therm2:
          ! we will use the Sig Wave Height from WW3, but could use this for testing without WW3
          ! if (tr_fsd) wave_sig_ht(i,j,iblk) = c4*SQRT(SUM(wave_spectrum(i,j,:,iblk)*dwavefreq(:)))

          call icepack_step_therm2(&
               dt=config_dt, &
               hin_max=categoryThicknessLimits(:), &
               aicen=iceAreaCategory(1,:,iCell), &
               vicen=iceVolumeCategory(1,:,iCell), &
               vsnon=snowVolumeCategory(1,:,iCell), &
               aicen_init=iceAreaCategoryInitial(:,iCell), &  
               vicen_init=iceVolumeCategoryInitial(:,iCell), &
               trcrn=tracerArrayCategory, &
               aice0=openWaterArea(iCell), &
               aice=iceAreaCell(iCell), &
               trcr_depend=ciceTracerObject % parentIndex,   &
               trcr_base=ciceTracerObject % firstAncestorMask, & 
               n_trcr_strata=ciceTracerObject % ancestorNumber, &
               nt_strata=ciceTracerObject % ancestorIndices, &
               Tf=seaFreezingTemperature(iCell), &
               sss=seaSurfaceSalinity(iCell), &
               salinz=initialSalinityProfile(:,iCell), &
               rsiden=lateralIceMeltFractionCategory(:,iCell),   &
               meltl=lateralIceMelt(iCell), &
               wlat=lateralIceMeltRate(iCell), &
               frzmlt=freezingMeltingPotential(iCell), &
               frazil=frazilFormation(iCell), &
               frain=rainfallRate(iCell), &
               fpond=pondFreshWaterFlux(iCell), &
               fresh=oceanFreshWaterFlux(iCell), &
               fsalt=oceanSaltFlux(iCell), &
               fhocn=oceanHeatFlux(iCell), &
               faero_ocn=oceanAerosolFlux(:,iCell), &
               first_ice=newlyFormedIceLogical(:), & 
               flux_bio=oceanBioFluxesTemp(:,iCell), &
               ocean_bio=oceanBioConcentrationsInUse(:,iCell), &
               frazil_diag=frazilGrowthDiagnostic(iCell), &
               frz_onset=freezeOnset(iCell), & ! optional
               yday=dayOfYear) ! optional
!               yday=dayOfYear, & ! optional
!               fiso_ocn=, & !optional - isotope flux to ocean  (kg/m^2/s)
!               HDO_ocn=, & ! optional - ocean concentration of HDO (kg/kg)
!               H2_16O_ocn=, & ! optional - ocean concentration of H2_16O (kg/kg) 
!               H2_18O_ocn=, & ! optional - ocean concentration of H2_18O (kg/kg)
!               nfsd=nFloeCategories, &
!               wave_sig_ht=SignificantWaveHeight(iCell), &
!               wave_spectrum=WaveSpectra(:,iCell),    &
!               wavefreq=WaveFrequency(:),       &
!               dwavefreq=WaveFreqBinWidth(:),      &
!               d_afsd_latg=DFloeSizeLateralGrowth(:,iCell),    &
!               d_afsd_newi=DFloeSizeNewIce(:,iCell),   &
!               d_afsd_latm=DFloeSizeLateralMelt(:,iCell),   &
!               d_afsd_weld=DFloeSizeWeld(:,iCell),   &
!               floe_rad_c=FloeSizeBinCenter(:), &
!               floe_binwidth=FloeSizeBinWidth(:))

          do iBioTracers = 1, ciceTracerObject % nBioTracers
             oceanBioFluxes(iBioTracers,iCell) = oceanBioFluxes(iBioTracers,iCell) + oceanBioFluxesTemp(iBioTracers,iCell)
          enddo

          ! update
          do iCategory = 1, nCategories
             newlyFormedIce(iCategory,iCell) = 0
             if (newlyFormedIceLogical(iCategory)) newlyFormedIce(iCategory,iCell) = 1
          enddo ! iCategory

          ! get category tracer array
          call get_cice_tracer_array_category(block, ciceTracerObject, &
                 tracerArrayCategory, iCell, setGetPhysicsTracers, setGetBGCTracers)

          abortFlag = icepack_warnings_aborted()
          ! cell-specific abort message
          if (abortFlag) then
             call mpas_log_write("column_itd_thermodynamics: "//trim(abortMessage) , messageType=MPAS_LOG_ERR)
             call mpas_log_write("iCell: $i", messageType=MPAS_LOG_ERR, intArgs=(/indexToCellID(iCell)/))
          endif
          call seaice_icepack_write_warnings(abortFlag)

       enddo ! iCell

       ! error checking
       call seaice_critical_error_write_block(domain, block, abortFlag)
       call seaice_check_critical_error(domain, abortFlag)

       ! newly formed ice
       deallocate(newlyFormedIceLogical)
       deallocate(oceanBioFluxesTemp)

       block => block % next
    end do

  end subroutine column_itd_thermodynamics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_prep_radiation
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine column_prep_radiation(domain)

    use icepack_intfc, only: &
         icepack_prep_radiation

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate, &
         atmos_coupling, &
         shortwave

    ! dimensions
    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nIceLayers, &
         nSnowLayers

    ! variables
    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown, &
         shortwaveScalingFactor, &
         albedoVisibleDirectArea, &
         albedoVisibleDiffuseArea, &
         albedoIRDirectArea, &
         albedoIRDiffuseArea

    real(kind=RKIND), dimension(:,:), pointer :: &
         surfaceShortwaveFlux, &
         interiorShortwaveFlux, &
         penetratingShortwaveFlux, &
         penetratingShortwaveFluxVisibleDirect, &
         penetratingShortwaveFluxVisibleDiffuse, &
         penetratingShortwaveFluxIRDirect, &
         penetratingShortwaveFluxIRDiffuse

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         shortwaveLayerPenetration, &
         absorbedShortwaveSnowLayer, &
         absorbedShortwaveIceLayer

    ! local
    integer :: &
         iCell

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)

       call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)

       call MPAS_pool_get_array(atmos_coupling, "shortwaveVisibleDirectDown", shortwaveVisibleDirectDown)
       call MPAS_pool_get_array(atmos_coupling, "shortwaveVisibleDiffuseDown", shortwaveVisibleDiffuseDown)
       call MPAS_pool_get_array(atmos_coupling, "shortwaveIRDirectDown", shortwaveIRDirectDown)
       call MPAS_pool_get_array(atmos_coupling, "shortwaveIRDiffuseDown", shortwaveIRDiffuseDown)

       call MPAS_pool_get_array(shortwave, "albedoVisibleDirectArea", albedoVisibleDirectArea)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDiffuseArea", albedoVisibleDiffuseArea)
       call MPAS_pool_get_array(shortwave, "albedoIRDirectArea", albedoIRDirectArea)
       call MPAS_pool_get_array(shortwave, "albedoIRDiffuseArea", albedoIRDiffuseArea)
       call MPAS_pool_get_array(shortwave, "shortwaveScalingFactor", shortwaveScalingFactor)
       call MPAS_pool_get_array(shortwave, "surfaceShortwaveFlux", surfaceShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "interiorShortwaveFlux", interiorShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "penetratingShortwaveFlux", penetratingShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "penetratingShortwaveFluxVisibleDirect", penetratingShortwaveFluxVisibleDirect)
       call MPAS_pool_get_array(shortwave, "penetratingShortwaveFluxVisibleDiffuse", penetratingShortwaveFluxVisibleDiffuse)
       call MPAS_pool_get_array(shortwave, "penetratingShortwaveFluxIRDirect", penetratingShortwaveFluxIRDirect)
       call MPAS_pool_get_array(shortwave, "penetratingShortwaveFluxIRDiffuse", penetratingShortwaveFluxIRDiffuse)
       call MPAS_pool_get_array(shortwave, "shortwaveLayerPenetration", shortwaveLayerPenetration)
       call MPAS_pool_get_array(shortwave, "absorbedShortwaveSnowLayer", absorbedShortwaveSnowLayer)
       call MPAS_pool_get_array(shortwave, "absorbedShortwaveIceLayer", absorbedShortwaveIceLayer)

       do iCell = 1, nCellsSolve

          call icepack_prep_radiation(&
               aice=iceAreaCell(iCell), &
               aicen=iceAreaCategory(1,:,iCell), &
               swvdr=shortwaveVisibleDirectDown(iCell), &
               swvdf=shortwaveVisibleDiffuseDown(iCell), &
               swidr=shortwaveIRDirectDown(iCell), &
               swidf=shortwaveIRDiffuseDown(iCell), &
               alvdr_ai=albedoVisibleDirectArea(iCell), &
               alvdf_ai=albedoVisibleDiffuseArea(iCell), &
               alidr_ai=albedoIRDirectArea(iCell), &
               alidf_ai=albedoIRDiffuseArea(iCell), &
               scale_factor=shortwaveScalingFactor(iCell), &
               fswsfcn=surfaceShortwaveFlux(:,iCell), &
               fswintn=interiorShortwaveFlux(:,iCell), &
               fswthrun=penetratingShortwaveFlux(:,iCell), &
               fswthrun_vdr=penetratingShortwaveFluxVisibleDirect(:,iCell), &
               fswthrun_vdf=penetratingShortwaveFluxVisibleDiffuse(:,iCell), &
               fswthrun_idr=penetratingShortwaveFluxIRDirect(:,iCell), &
               fswthrun_idf=penetratingShortwaveFluxIRDiffuse(:,iCell), &
               fswpenln=shortwaveLayerPenetration(:,:,iCell), &
               Sswabsn=absorbedShortwaveSnowLayer(:,:,iCell), &
               Iswabsn=absorbedShortwaveIceLayer(:,:,iCell))

       enddo ! iCell

       block => block % next
    end do

    call seaice_icepack_write_warnings(icepack_warnings_aborted())

  end subroutine column_prep_radiation

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_snow
!
!> \brief Enable snow grain aging, effective snow density, wind compaction and redistribution
!>
!> \author Nicole Jeffery, LANL
!> \date 3rd April 2017. Modified for icepack July 2023
!> \details
!>
!> Snow physics improvements include:
!> 1) Snow redistribution by wind (multiple options available).
!> Includes parametrizations for snow compaction, redistribution among categories/level ice
!> and loss to leads.
!> 2) Tracking of snow liquid and ice content.  Liquid is stored in snow before passing to ponds.
!> Effective snow density is also tracked.
!> 3) Snow grain radius aging based on wet (liquid content) and dry (temperature gradient) metamorphism.
!> 4) Effective snow density (based on snow liquid/ice content and compaction)
!
!-----------------------------------------------------------------------

  subroutine column_snow(domain)

    use icepack_intfc, only: &
         icepack_step_snow, &
         icepack_warnings_clear, &
         icepack_warnings_aborted

    use seaice_constants, only: &
         seaicePuny, &
         seaiceDensitySnow

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate, &
         atmos_forcing, &
         snow, &
         ocean_fluxes, &
         atmos_coupling

    logical, pointer :: &
         config_use_effective_snow_density, &
         config_use_snow_grain_radius, &
         config_use_column_biogeochemistry

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         snowIceMass, &
         snowLiquidMass, &
         snowDensity, &
         snowVolumeCategory, &
         iceAreaCategory, &
         iceVolumeCategory, &
         levelIceArea, &
         levelIceVolume, &
         iceEnthalpy, &
         snowEnthalpy, &
         iceSalinity, &
         surfaceTemperature, &
         snowGrainRadius, &
         snowEmpiricalGrowthParameterTau, &
         snowEmpiricalGrowthParameterKappa, &
         snowPropertyRate

    real(kind=RKIND), dimension(:,:), pointer :: &
         snowMeltMassCategory

    real(kind=RKIND), dimension(:), pointer :: &
         windSpeed, &
         oceanFreshWaterFlux, &
         oceanHeatFlux, &
         snowLossToLeads, &
         snowfallRate, &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         snowDensityViaContent, &
         snowDensityViaCompaction, &
         snowMeltMassCell

    real(kind=RKIND), pointer :: &
         config_dt

    integer, pointer :: &
         nCellsSolve, &
         nSnowLayers, &
         nIceLayers, &
         nCategories

    integer :: &
         iCell, &
         iSnowLayer, &
         iIceLayer, &
         iCategory

    logical :: &
         abortFlag, &
         setGetPhysicsTracers, &
         setGetBGCTracers

    character(len=strKIND) :: &
         abortMessage, &
         abortLocation

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
       call MPAS_pool_get_subpool(block % structs, "snow", snow)
       call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmos_forcing)
       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", ocean_fluxes)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)

       call MPAS_pool_get_config(block % configs, "config_use_effective_snow_density", config_use_effective_snow_density)
       call MPAS_pool_get_config(block % configs, "config_use_snow_grain_radius", config_use_snow_grain_radius)
       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nCategories", nCategories)
       call MPAS_pool_get_dimension(block % dimensions, "nSnowLayers", nSnowLayers)
       call MPAS_pool_get_dimension(block % dimensions, "nIceLayers", nIceLayers)

       call MPAS_pool_get_array(snow, "snowDensityViaContent", snowDensityViaContent)
       call MPAS_pool_get_array(snow, "snowDensityViaCompaction", snowDensityViaCompaction)
       call MPAS_pool_get_array(snow, "snowMeltMassCategory", snowMeltMassCategory)
       call MPAS_pool_get_array(snow, "snowMeltMassCell", snowMeltMassCell)
       call MPAS_pool_get_array(snow, "snowLossToLeads", snowLossToLeads)

       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracers, "snowIceMass", snowIceMass, 1)
       call MPAS_pool_get_array(tracers, "snowLiquidMass", snowLiquidMass, 1)
       call MPAS_pool_get_array(tracers, "snowDensity", snowDensity, 1)
       call MPAS_pool_get_array(tracers, "snowGrainRadius", snowGrainRadius, 1)
       call MPAS_pool_get_array(tracers, "levelIceArea", levelIceArea, 1)
       call MPAS_pool_get_array(tracers, "levelIceVolume", levelIceVolume, 1)
       call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy, 1)
       call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy, 1)
       call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity, 1)
       call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

       call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(tracers_aggregate, "snowVolumeCell", snowVolumeCell)

       call MPAS_pool_get_array(atmos_coupling, "snowfallRate", snowfallRate)

       call MPAS_pool_get_array(atmos_forcing, "windSpeed", windSpeed)

       call MPAS_pool_get_array(ocean_fluxes, "oceanFreshWaterFlux", oceanFreshWaterFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanHeatFlux", oceanHeatFlux)

       setGetPhysicsTracers = .true.
       setGetBGCTracers     = config_use_column_biogeochemistry

       do iCell = 1, nCellsSolve

          abortFlag = .false.
          abortMessage = ""

          call icepack_warnings_clear()
          call icepack_step_snow (&
               dt=config_dt, &
               wind=windSpeed(iCell), &
               aice=iceAreaCell(iCell), &
               aicen=iceAreaCategory(1,:,iCell), &
               vicen=iceVolumeCategory(1,:,iCell), &
               vsnon=snowVolumeCategory(1,:,iCell), &
               alvl=levelIceArea(1,:,iCell), &
               vlvl=levelIceVolume(1,:,iCell), &
               smice=snowIceMass(:,:,iCell), &
               smliq=snowLiquidMass(:,:,iCell), &
               rhos_cmpn=snowDensity(:,:,iCell), &
               rsnw=snowGrainRadius(:,:,iCell), &
               zqin1=iceEnthalpy(1,:,iCell), &
               zSin1=iceSalinity(1,:,iCell), &
               Tsfc=surfaceTemperature(1,:,iCell), &
               zqsn=snowEnthalpy(:,:,iCell), &
               fresh=oceanFreshWaterFlux(iCell), &
               fhocn=oceanHeatFlux(iCell), &
               fsloss=snowLossToLeads(iCell), &
               fsnow=snowfallRate(iCell) &
               )

          abortFlag = icepack_warnings_aborted()
          call seaice_icepack_write_warnings(abortFlag)

          if (config_use_snow_grain_radius) then
             snowDensityViaContent(iCell) = 0.0_RKIND
             do iCategory = 1, nCategories
                if (snowVolumeCategory(1,iCategory,iCell) .gt. 0.0_RKIND) then
                   do iSnowLayer = 1, nSnowLayers
                      snowDensityViaContent(iCell) = snowDensityViaContent(iCell) &
                         + snowVolumeCategory(1,iCategory,iCell) * &
                         (snowIceMass(iSnowLayer,iCategory,iCell) + &
                         snowLiquidMass(iSnowLayer,iCategory,iCell))
                   enddo !iSnowLayer
                endif !snowVolumeCategory
             enddo !iCategory
             if (snowVolumeCell(iCell) .gt.  seaicePuny) then
                snowDensityViaContent(iCell) = snowDensityViaContent(iCell)/ &
                   (snowVolumeCell(iCell) * real(nSnowLayers,kind=RKIND))  !!!CHECK THIS!!!
             else
                snowDensityViaContent(iCell)    = seaiceDensitySnow
             endif !snowVolumeCell
          endif

       enddo !iCell

       block => block % next
    end do

  end subroutine column_snow

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_radiation
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine column_radiation(domain, clock, lInitialization)

    use icepack_intfc, only: &
         icepack_step_radiation

    use seaice_constants, only: &
         seaicePi

    type(domain_type), intent(inout) :: domain

    type(MPAS_clock_type), intent(in) :: clock

    logical, intent(in) :: &
         lInitialization

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         atmos_coupling, &
         shortwave, &
         ponds, &
         biogeochemistry

    ! configs
    real(kind=RKIND), pointer :: &
         config_dt

    logical, pointer :: &
         config_use_snow_grain_radius, &
         config_use_shortwave_bioabsorption, &
         config_use_brine, &
         config_use_modal_aerosols, &
         config_use_column_biogeochemistry, &
         config_use_zaerosols

    ! dimensions
    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nIceLayers, &
         nSnowLayers, &
         nAerosols, &
         nAlgae, &
         nBioLayers, &
         nzAerosols, &
         maxAerosolType

    ! variables
    real(kind=RKIND), dimension(:), pointer :: &
         latCell, &
         lonCell, &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown, &
         solarZenithAngleCosine, &
         snowfallRate

    real(kind=RKIND), dimension(:,:), pointer :: &
         surfaceShortwaveFlux, &
         interiorShortwaveFlux, &
         penetratingShortwaveFlux, &
         penetratingShortwaveFluxVisibleDirect, &
         penetratingShortwaveFluxVisibleDiffuse, &
         penetratingShortwaveFluxIRDirect, &
         penetratingShortwaveFluxIRDiffuse, &
         bareIceAlbedoCategory, &
         snowAlbedoCategory, &
         pondAlbedoCategory, &
         effectivePondAreaCategory, &
         pondSnowDepthDifference, &
         pondLidMeltFluxFraction, &
         albedoVisibleDirectCategory, &
         albedoVisibleDiffuseCategory, &
         albedoIRDirectCategory, &
         albedoIRDiffuseCategory, &
         snowFractionCategory

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature, &
         levelIceArea, &
         pondArea, &
         pondDepth, &
         pondLidThickness, &
         shortwaveLayerPenetration, &
         absorbedShortwaveSnowLayer, &
         absorbedShortwaveIceLayer, &
         snowScatteringAerosol, &
         snowBodyAerosol, &
         iceScatteringAerosol, &
         iceBodyAerosol, &
         brineFraction, &
         bioTracerShortwave, &
         snowGrainRadius, &
         verticalAerosolsConc, &
         verticalAlgaeConc

    real(kind=RKIND), pointer :: &
         dayOfNextShortwaveCalculation ! needed for CESM like coupled simulations

    character(len=strKIND), pointer :: &
         config_calendar_type

    character(len=strKIND) :: &
         calendarType ! needed for CESM like coupled simulations

    ! aerosols and snow grain radius arrays
    real(kind=RKIND), dimension(:,:), allocatable :: &
         aerosolsArray, &
         snow_grain_radius

    ! local
    integer :: &
         iCell, &
         iSnowLayer, &
         iCategory, &
         iAerosol, &
         iBioTracers

    real(kind=RKIND) :: &
         dayOfYear, &
         lonCellColumn

    integer :: &
         secondsIntoDay, &
         daysInYear

    logical :: &
         setGetPhysicsTracers, &
         setGetBGCTracers

    ! day of year
    call get_day_of_year(clock, dayOfYear)

    ! seconds into day
    call get_seconds_into_day(clock, secondsIntoDay)

    ! get days in year
    call get_days_in_year(domain, clock, daysInYear)

    call MPAS_pool_get_config(domain % configs, "config_use_snow_grain_radius", config_use_snow_grain_radius)
    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)
    call MPAS_pool_get_config(domain % configs, "config_use_shortwave_bioabsorption", config_use_shortwave_bioabsorption)
    call MPAS_pool_get_config(domain % configs, "config_use_modal_aerosols",config_use_modal_aerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry",config_use_column_biogeochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_zaerosols",config_use_zaerosols)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "ponds", ponds)
       call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nCategories", nCategories)
       call MPAS_pool_get_dimension(block % dimensions, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(block % dimensions, "nSnowLayers", nSnowLayers)
       call MPAS_pool_get_dimension(block % dimensions, "nAerosols", nAerosols)
       call MPAS_pool_get_dimension(block % dimensions, "nAlgae", nAlgae)
       call MPAS_pool_get_dimension(block % dimensions, "nBioLayers", nBioLayers)
       call MPAS_pool_get_dimension(block % dimensions, "nzAerosols", nzAerosols)
       call MPAS_pool_get_dimension(block % dimensions, "maxAerosolType", maxAerosolType)

       call MPAS_pool_get_array(mesh, "latCell", latCell)
       call MPAS_pool_get_array(mesh, "lonCell", lonCell)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
       call MPAS_pool_get_array(tracers, "levelIceArea", levelIceArea, 1)
       call MPAS_pool_get_array(tracers, "pondArea", pondArea, 1)
       call MPAS_pool_get_array(tracers, "pondDepth", pondDepth, 1)
       call MPAS_pool_get_array(tracers, "pondLidThickness", pondLidThickness, 1)
       call MPAS_pool_get_array(tracers, "snowScatteringAerosol", snowScatteringAerosol, 1)
       call MPAS_pool_get_array(tracers, "snowBodyAerosol", snowBodyAerosol, 1)
       call MPAS_pool_get_array(tracers, "iceScatteringAerosol", iceScatteringAerosol, 1)
       call MPAS_pool_get_array(tracers, "iceBodyAerosol", iceBodyAerosol, 1)
       call MPAS_pool_get_array(tracers, "brineFraction", brineFraction, 1)
       call MPAS_pool_get_array(tracers, "snowGrainRadius", snowGrainRadius, 1)
       call MPAS_pool_get_array(tracers, "verticalAlgaeConc", verticalAlgaeConc, 1)
       call MPAS_pool_get_array(tracers, "verticalAerosolsConc", verticalAerosolsConc, 1)

       call MPAS_pool_get_array(atmos_coupling, "shortwaveVisibleDirectDown", shortwaveVisibleDirectDown)
       call MPAS_pool_get_array(atmos_coupling, "shortwaveVisibleDiffuseDown", shortwaveVisibleDiffuseDown)
       call MPAS_pool_get_array(atmos_coupling, "shortwaveIRDirectDown", shortwaveIRDirectDown)
       call MPAS_pool_get_array(atmos_coupling, "shortwaveIRDiffuseDown", shortwaveIRDiffuseDown)
       call MPAS_pool_get_array(atmos_coupling, "snowfallRate", snowfallRate)

       call MPAS_pool_get_array(shortwave, "dayOfNextShortwaveCalculation", dayOfNextShortwaveCalculation)
       call MPAS_pool_get_array(shortwave, "solarZenithAngleCosine", solarZenithAngleCosine)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDirectCategory", albedoVisibleDirectCategory)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDiffuseCategory", albedoVisibleDiffuseCategory)
       call MPAS_pool_get_array(shortwave, "albedoIRDirectCategory", albedoIRDirectCategory)
       call MPAS_pool_get_array(shortwave, "albedoIRDiffuseCategory", albedoIRDiffuseCategory)
       call MPAS_pool_get_array(shortwave, "surfaceShortwaveFlux", surfaceShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "interiorShortwaveFlux", interiorShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "penetratingShortwaveFlux", penetratingShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "penetratingShortwaveFluxVisibleDirect", penetratingShortwaveFluxVisibleDirect)
       call MPAS_pool_get_array(shortwave, "penetratingShortwaveFluxVisibleDiffuse", penetratingShortwaveFluxVisibleDiffuse)
       call MPAS_pool_get_array(shortwave, "penetratingShortwaveFluxIRDirect", penetratingShortwaveFluxIRDirect)
       call MPAS_pool_get_array(shortwave, "penetratingShortwaveFluxIRDiffuse", penetratingShortwaveFluxIRDiffuse)
       call MPAS_pool_get_array(shortwave, "shortwaveLayerPenetration", shortwaveLayerPenetration)
       call MPAS_pool_get_array(shortwave, "absorbedShortwaveSnowLayer", absorbedShortwaveSnowLayer)
       call MPAS_pool_get_array(shortwave, "absorbedShortwaveIceLayer", absorbedShortwaveIceLayer)
       call MPAS_pool_get_array(shortwave, "bareIceAlbedoCategory", bareIceAlbedoCategory)
       call MPAS_pool_get_array(shortwave, "snowAlbedoCategory", snowAlbedoCategory)
       call MPAS_pool_get_array(shortwave, "pondAlbedoCategory", pondAlbedoCategory)
       call MPAS_pool_get_array(shortwave, "effectivePondAreaCategory", effectivePondAreaCategory)
       call MPAS_pool_get_array(shortwave, "snowFractionCategory", snowFractionCategory)

       call MPAS_pool_get_array(ponds, "pondSnowDepthDifference", pondSnowDepthDifference)
       call MPAS_pool_get_array(ponds, "pondLidMeltFluxFraction", pondLidMeltFluxFraction)

       call MPAS_pool_get_array(biogeochemistry, "bioTracerShortwave", bioTracerShortwave)

       ! calendar type
       call MPAS_pool_get_config(block % configs, "config_calendar_type", config_calendar_type)
       if (trim(config_calendar_type) == "gregorian") then
          calendarType = "GREGORIAN"
       else
          calendarType = "NOLEAP"
       endif

       ! aerosols array
       allocate(aerosolsArray(4*nAerosols,nCategories))

       setGetPhysicsTracers = .true.
       setGetBGCTracers     = (config_use_column_biogeochemistry .or. config_use_zaerosols)

       ! snow grain radius array
       allocate(snow_grain_radius(nSnowLayers,nCategories))

       !$omp parallel do default(shared) &
       !$omp&   firstprivate(aerosolsArray,snow_grain_radius) &
       !$omp&   private(iCategory,iAerosol,lonCellColumn,iSnowLayer)

       do iCell = 1, nCellsSolve

          ! set aerosols array
          do iCategory = 1, nCategories
             do iAerosol = 1, nAerosols

                aerosolsArray(1+4*(iAerosol-1), iCategory) = snowScatteringAerosol(iAerosol,iCategory,iCell)
                aerosolsArray(2+4*(iAerosol-1), iCategory) = snowBodyAerosol(iAerosol,iCategory,iCell)
                aerosolsArray(3+4*(iAerosol-1), iCategory) = iceScatteringAerosol(iAerosol,iCategory,iCell)
                aerosolsArray(4+4*(iAerosol-1), iCategory) = iceBodyAerosol(iAerosol,iCategory,iCell)

             enddo ! iAerosol
          enddo ! iCategory

          snow_grain_radius(:,:) = 0.0_RKIND
          if (config_use_snow_grain_radius) then
             do iCategory = 1, nCategories
                do iSnowLayer = 1, nSnowLayers
                   snow_grain_radius(iSnowLayer, iCategory) = snowGrainRadius(iSnowLayer,iCategory,iCell)
                enddo ! iSnowLayer
             enddo ! iCategory
          endif

          lonCellColumn = lonCell(iCell)
          if (lonCellColumn > seaicePi) lonCellColumn = lonCellColumn - 2.0_RKIND * seaicePi

          call icepack_step_radiation(&
               dt=config_dt, &
               fbri=brineFraction(1,:,iCell), &
               aicen=iceAreaCategory(1,:,iCell), &
               vicen=iceVolumeCategory(1,:,iCell), &
               vsnon=snowVolumeCategory(1,:,iCell), &
               Tsfcn=surfaceTemperature(1,:,iCell), &
               alvln=levelIceArea(1,:,iCell), &
               apndn=pondArea(1,:,iCell), &
               hpndn=pondDepth(1,:,iCell), &
               ipndn=pondLidThickness(1,:,iCell), &
               aeron=aerosolsArray, &
               bgcNn=verticalAlgaeConc(:,:,iCell), &
               zaeron=verticalAerosolsConc(:,:,iCell), &
               trcrn_bgcsw=bioTracerShortwave(:,:,iCell), &
               TLAT=latCell(iCell), &
               TLON=lonCellColumn, &
               calendar_type=calendarType, &
               days_per_year=daysInYear, &
               nextsw_cday=dayOfNextShortwaveCalculation, &
               yday=dayOfYear, &
               sec=secondsIntoDay, &
               swvdr=shortwaveVisibleDirectDown(iCell), &
               swvdf=shortwaveVisibleDiffuseDown(iCell), &
               swidr=shortwaveIRDirectDown(iCell), &
               swidf=shortwaveIRDiffuseDown(iCell), &
               coszen=solarZenithAngleCosine(iCell), &
               fsnow=snowfallRate(iCell), &
               alvdrn=albedoVisibleDirectCategory(:,iCell), &
               alvdfn=albedoVisibleDiffuseCategory(:,iCell), &
               alidrn=albedoIRDirectCategory(:,iCell), &
               alidfn=albedoIRDiffuseCategory(:,iCell), &
               fswsfcn=surfaceShortwaveFlux(:,iCell), &
               fswintn=interiorShortwaveFlux(:,iCell), &
               fswthrun=penetratingShortwaveFlux(:,iCell), &
               fswthrun_vdr=penetratingShortwaveFluxVisibleDirect(:,iCell), &
               fswthrun_vdf=penetratingShortwaveFluxVisibleDiffuse(:,iCell), &
               fswthrun_idr=penetratingShortwaveFluxIRDirect(:,iCell), &
               fswthrun_idf=penetratingShortwaveFluxIRDiffuse(:,iCell), &
               fswpenln=shortwaveLayerPenetration(:,:,iCell), &
               Sswabsn=absorbedShortwaveSnowLayer(:,:,iCell), &
               Iswabsn=absorbedShortwaveIceLayer(:,:,iCell), &
               albicen=bareIceAlbedoCategory(:,iCell), &
               albsnon=snowAlbedoCategory(:,iCell), &
               albpndn=pondAlbedoCategory(:,iCell), &
               apeffn=effectivePondAreaCategory(:,iCell), &
               snowfracn=snowFractionCategory(:,iCell), &
               dhsn=pondSnowDepthDifference(:,iCell), &
               ffracn=pondLidMeltFluxFraction(:,iCell), &
               rsnow=snow_grain_radius(:,:), &
               l_print_point=.false., &
               initonly=lInitialization)

       enddo ! iCell

       deallocate(snow_grain_radius)

       ! aerosols array
       deallocate(aerosolsArray)

       block => block % next
    end do

    call seaice_icepack_write_warnings(icepack_warnings_aborted())

  end subroutine column_radiation

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_ridging
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine column_ridging(domain)

    use icepack_intfc, only: &
         icepack_step_ridge

    use seaice_constants, only: &
         seaicePuny

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         icestate, &
         tracers, &
         tracers_aggregate, &
         ponds, &
         ocean_fluxes, &
         ocean_coupling, &
         ridging, &
         aerosols, &
         biogeochemistry, &
         initial, &
         velocity_solver

    ! configs
    logical, pointer :: &
         config_use_column_biogeochemistry, &
         config_use_zaerosols

    real(kind=RKIND), pointer :: &
         config_dt

    integer, pointer :: &
         config_dynamics_subcycle_number

    ! dimensions
    integer, pointer :: &
         nCellsSolve, &
         nCategories

    ! variables
    real(kind=RKIND), dimension(:), pointer :: &
         pondFreshWaterFlux, &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         seaFreezingTemperature, &
         iceAreaCell, &
         ridgeConvergence, &
         ridgeShear, &
         openWaterArea, &
         areaLossRidge, &
         areaGainRidge, &
         iceVolumeRidged, &
         openingRateRidge, &
         categoryThicknessLimits

    real(kind=RKIND), dimension(:,:), pointer :: &
         oceanAerosolFlux, &
         ridgeParticipationFunction, &
         ratioRidgeThicknessToIce, &
         fractionNewRidgeArea, &
         fractionNewRidgeVolume, &
         areaLossRidgeCategory, &
         areaGainRidgeCategory, &
         iceVolumeRidgedCategory, &
         raftingIceArea, &
         raftingIceVolume, &
         oceanBioFluxes

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    integer, dimension(:,:), pointer :: &
         newlyFormedIce

    integer, dimension(:), pointer :: &
         indexToCellID

    real(kind=RKIND), pointer :: &
         dynamicsTimeStep

    ! local
    integer :: &
         iCell, &
         iCategory, &
         iBioTracers, &
         iBioData, &
         iBioLayers

    real(kind=RKIND), dimension(:,:), allocatable :: &
         oceanBioFluxesTemp

    logical, dimension(:), allocatable :: &
         newlyFormedIceLogical

    logical :: &
         setGetPhysicsTracers, &
         setGetBGCTracers

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
       call MPAS_pool_get_subpool(block % structs, "ponds", ponds)
       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", ocean_fluxes)
       call MPAS_pool_get_subpool(block % structs, "ridging", ridging)
       call MPAS_pool_get_subpool(block % structs, "aerosols", aerosols)
       call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)
       call MPAS_pool_get_subpool(block % structs, "initial", initial)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

       call MPAS_pool_get_config(block % configs, "config_dynamics_subcycle_number", config_dynamics_subcycle_number)
       call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
       call MPAS_pool_get_config(block % configs, "config_use_zaerosols", config_use_zaerosols)
       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)

       call MPAS_pool_get_array(velocity_solver, "dynamicsTimeStep", dynamicsTimeStep)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

       call MPAS_pool_get_array(mesh, "indexToCellID", indexToCellID)

       call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)

       call MPAS_pool_get_array(icestate, "openWaterArea", openWaterArea)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

       call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)

       call MPAS_pool_get_array(ocean_fluxes, "oceanFreshWaterFlux", oceanFreshWaterFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanSaltFlux", oceanSaltFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanHeatFlux", oceanHeatFlux)

       call MPAS_pool_get_array(ridging, "ridgeConvergence", ridgeConvergence)
       call MPAS_pool_get_array(ridging, "ridgeShear", ridgeShear)
       call MPAS_pool_get_array(ridging, "areaLossRidge", areaLossRidge)
       call MPAS_pool_get_array(ridging, "areaGainRidge", areaGainRidge)
       call MPAS_pool_get_array(ridging, "iceVolumeRidged", iceVolumeRidged)
       call MPAS_pool_get_array(ridging, "openingRateRidge", openingRateRidge)
       call MPAS_pool_get_array(ridging, "ridgeParticipationFunction", ridgeParticipationFunction)
       call MPAS_pool_get_array(ridging, "ratioRidgeThicknessToIce", ratioRidgeThicknessToIce)
       call MPAS_pool_get_array(ridging, "fractionNewRidgeArea", fractionNewRidgeArea)
       call MPAS_pool_get_array(ridging, "fractionNewRidgeVolume", fractionNewRidgeVolume)
       call MPAS_pool_get_array(ridging, "areaLossRidgeCategory", areaLossRidgeCategory)
       call MPAS_pool_get_array(ridging, "areaGainRidgeCategory", areaGainRidgeCategory)
       call MPAS_pool_get_array(ridging, "iceVolumeRidgedCategory", iceVolumeRidgedCategory)
       call MPAS_pool_get_array(ridging, "raftingIceArea", raftingIceArea)
       call MPAS_pool_get_array(ridging, "raftingIceVolume", raftingIceVolume)

       call MPAS_pool_get_array(aerosols, "oceanAerosolFlux", oceanAerosolFlux)

       call MPAS_pool_get_array(ponds, "pondFreshWaterFlux", pondFreshWaterFlux)

       call MPAS_pool_get_array(biogeochemistry, "newlyFormedIce", newlyFormedIce)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioFluxes", oceanBioFluxes)

       call MPAS_pool_get_array(initial, "categoryThicknessLimits", categoryThicknessLimits)

       ! newly formed ice
       allocate(newlyFormedIceLogical(nCategories))

       setGetPhysicsTracers = .true.
       setGetBGCTracers     = (config_use_column_biogeochemistry .or. config_use_zaerosols)

       allocate(oceanBioFluxesTemp(ciceTracerObject % nBioTracers,nCellsSolve))

       do iCell = 1, nCellsSolve

          ! newly formed ice
          do iCategory = 1, nCategories
             newlyFormedIceLogical(iCategory) = (newlyFormedIce(iCategory,iCell) == 1)
          enddo ! iCategory

          ! set the category tracer array
          call set_cice_tracer_array_category(block, ciceTracerObject, &
                 tracerArrayCategory, iCell, setGetPhysicsTracers, setGetBGCTracers)

          oceanBioFluxesTemp(:,iCell) = 0.0_RKIND

          call icepack_step_ridge(&
               dt=dynamicsTimeStep, &
               ndtd=config_dynamics_subcycle_number, &
               hin_max=categoryThicknessLimits, & ! hin_max, dimension(0:ncat), intent(inout)
               rdg_conv=ridgeConvergence(iCell), &
               rdg_shear=ridgeShear(iCell), &
               aicen=iceAreaCategory(1,:,iCell), &
               trcrn=tracerArrayCategory, &              ! trcrn, dimension(:,:), intent(inout)
               vicen=iceVolumeCategory(1,:,iCell), &
               vsnon=snowVolumeCategory(1,:,iCell), &
               aice0=openWaterArea(iCell), &
               trcr_depend=ciceTracerObject % parentIndex, & ! trcr_depend
               trcr_base=ciceTracerObject % firstAncestorMask, & ! trcr_base
               n_trcr_strata=ciceTracerObject % ancestorNumber, & ! n_trcr_strata
               nt_strata=ciceTracerObject % ancestorIndices, & ! nt_strata
               dardg1dt=areaLossRidge(iCell), &
               dardg2dt=areaGainRidge(iCell), &
               dvirdgdt=iceVolumeRidged(iCell), &
               opening=openingRateRidge(iCell), &
               fpond=pondFreshWaterFlux(iCell), &
               fresh=oceanFreshWaterFlux(iCell), &
               fhocn=oceanHeatFlux(iCell), &
               faero_ocn=oceanAerosolFlux(:,iCell), & ! DC no fiso_ocn argument
               aparticn=ridgeParticipationFunction(:,iCell), &
               krdgn=ratioRidgeThicknessToIce(:,iCell), &
               aredistn=fractionNewRidgeArea(:,iCell), &
               vredistn=fractionNewRidgeVolume(:,iCell), &
               dardg1ndt=areaLossRidgeCategory(:,iCell), &
               dardg2ndt=areaGainRidgeCategory(:,iCell), &
               dvirdgndt=iceVolumeRidgedCategory(:,iCell), &
               araftn=raftingIceArea(:,iCell), &
               vraftn=raftingIceVolume(:,iCell), &
               aice=iceAreaCell(iCell), &
               fsalt=oceanSaltFlux(iCell), &
               first_ice=newlyFormedIceLogical(:), &
               flux_bio=oceanBioFluxesTemp(:,iCell), &
               Tf=seaFreezingTemperature(iCell))

          ! update
          do iCategory = 1, nCategories
             newlyFormedIce(iCategory,iCell) = 0
             if (newlyFormedIceLogical(iCategory)) newlyFormedIce(iCategory,iCell) = 1
          enddo ! iCategory

          do iBioTracers = 1, ciceTracerObject % nBioTracers
             oceanBioFluxes(iBioTracers,iCell) = oceanBioFluxes(iBioTracers,iCell) + oceanBioFluxesTemp(iBioTracers,iCell)
          enddo

          ! get category tracer array
          call get_cice_tracer_array_category(block, ciceTracerObject, &
                 tracerArrayCategory, iCell, setGetPhysicsTracers, setGetBGCTracers)

          call seaice_icepack_write_warnings(icepack_warnings_aborted())
       enddo ! iCell

       ! newly formed ice
       deallocate(newlyFormedIceLogical)
       deallocate(oceanBioFluxesTemp)

       block => block % next
    end do

  end subroutine column_ridging

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_biogeochemistry
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 19th October 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine column_biogeochemistry(domain)

    use icepack_intfc, only: &
         icepack_biogeochemistry, &
         icepack_load_ocean_bio_array, &
         icepack_warnings_clear

    use seaice_constants, only: &
         seaicePuny, &
         microgramsPerKilograms, &
         gramsIronPerMolIron

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         biogeochemistry, &
         diagnostics_biogeochemistry, &
         icestate, &
         tracers, &
         shortwave, &
         melt_growth_rates, &
         ocean_coupling, &
         atmos_coupling, &
         atmos_forcing, &
         initial

    ! configs
    real(kind=RKIND), pointer :: &
         config_dt

    logical, pointer :: &
         config_use_brine, &
         config_use_skeletal_biochemistry, &
         config_use_column_biogeochemistry, &
         config_use_zaerosols, &
         config_use_vertical_tracers, &
         config_use_atm_dust_file, &
         config_use_iron_solubility_file

    ! dimensions
    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nzAerosols, &
         nBioLayersP1, &
         nAlgae, &
         maxBCType, &
         maxDustType, &
         maxAerosolType

    ! variables

    real(kind=RKIND), dimension(:), pointer :: &
         netNitrateUptake, &
         netAmmoniumUptake, &
         netSpecificAlgalGrowthRate, &
         primaryProduction, &
         netBrineHeight, &
         seaSurfaceTemperature, &
         seaSurfaceSalinity, &
         seaFreezingTemperature, &
         snowfallRate, &
         totalSkeletalAlgae, &
         oceanNitrateConc, &
         oceanSilicateConc, &
         oceanAmmoniumConc, &
         oceanDMSConc, &
         oceanDMSPConc, &
         oceanHumicsConc, &
         openWaterArea, &
         totalChlorophyll, &
         IRON_Zolubility_wet, &
         IRON_Zolubility_dry, &
         IRON_in_duzt_fraction

    real(kind=RKIND), dimension(:,:), pointer :: &
         iceAreaCategoryInitial, &
         iceVolumeCategoryInitial, &
         snowVolumeCategoryInitial, &
         iceThicknessCategoryInitial, &  !icestate
         brineBottomChange, &
         brineTopChange, &
         darcyVelocityBio, &
         snowIceBioFluxes, &
         atmosIceBioFluxes, &
         oceanBioConcentrations, &
         oceanBioConcentrationsInUse, &
         oceanBioToIceInUse, &
         totalVerticalBiologyIce, &
         totalVerticalBiologySnow, &
         penetratingShortwaveFlux, &
         basalIceMeltCategory, &
         surfaceIceMeltCategory, &
         congelationCategory, &
         snowiceFormationCategory, &
         snowMeltCategory, &
         initialSalinityProfile, &
         atmosBioFluxes, &
         atmosBlackCarbonFlux, &
         atmosDustFlux, &
         atmosWetDustFlux, &
         atmosDryDustFlux, &
         oceanBioFluxes, &
         oceanAlgaeConc, &
         oceanDOCConc, &
         oceanDICConc, &
         oceanDONConc, &
         oceanParticulateIronConc, &
         oceanDissolvedIronConc, &
         oceanZAerosolConc, &
         bioShortwaveFluxCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         shortwaveLayerPenetration, &
         verticalNitrogenLosses, &
         bioPorosity, &
         bioTemperature, &
         bioDiffusivity, &
         bioPermeability, &
         bioShortwaveFlux, &
         iceAreaCategory, &    ! tracers  (1,ncat,ncell)
         iceVolumeCategory, &  ! tracers  (1,ncat,ncell)
         snowVolumeCategory, & ! tracers  (1,ncat,ncell)
         skeletalAlgaeConc, &
         oceanBioFluxesCategory, &
         brineFraction

    real(kind=RKIND), dimension(:,:), pointer :: &
         bgridTemperatureIceCell, &
         bgridSalinityIceCell, &
         bgridPorosityIceCell

    integer, dimension(:,:), pointer :: &
         newlyFormedIce

    integer, dimension(:), pointer :: &
         indexToCellID

    ! local
    integer :: &
         iCell, &
         iTracers, &
         iBioTracers, &
         iCategory, &
         iAlgae, &
         iBioData, &
         iBioCount, &
         iSnowCount, &
         iIceCount, &
         indexj, &
         iBioLayers, &
         iWarning, &
         nWarnings

    ! test carbon conservation
    real(kind=RKIND), dimension(:,:), allocatable :: &
         totalCarbonCatFinal, &
         totalCarbonCatInitial, &
         totalCarbonCatFlux, &
         brineHeightCatInitial, &
         brineHeightCatFinal

    real(kind=RKIND), dimension(:), allocatable :: &
         totalCarbonFinal, &
         totalCarbonInitial, &
         totalCarbonFlux

    real(kind=RKIND) :: &
         errorCheck

    logical, dimension(:), allocatable :: &
         newlyFormedIceLogical

    logical :: &
         abortFlag, &
         setGetPhysicsTracers, &
         setGetBGCTracers, &
         checkCarbon

    character(len=strKIND) :: &
         abortMessage, &
         abortLocation

    real(kind=RKIND), dimension(:), allocatable :: &
         carbonError

    real(kind=RKIND), parameter :: &
         accuracy = 1.0e-13_RKIND

    ! test carbon conservation
    checkCarbon = .false.

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)
       call MPAS_pool_get_subpool(block % structs, "diagnostics_biogeochemistry", diagnostics_biogeochemistry)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
       call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmos_forcing)
       call MPAS_pool_get_subpool(block % structs, "melt_growth_rates", melt_growth_rates)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
       call MPAS_pool_get_subpool(block % structs, "initial", initial)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nzAerosols", nzAerosols)
       call MPAS_pool_get_dimension(mesh, "nBioLayersP1", nBioLayersP1)
       call MPAS_pool_get_dimension(mesh, "nAlgae", nAlgae)
       call MPAS_pool_get_dimension(mesh, "maxBCType", maxBCType)
       call MPAS_pool_get_dimension(mesh, "maxDustType", maxDustType)

       call MPAS_pool_get_array(mesh, "indexToCellID", indexToCellID)

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_use_brine", config_use_brine)
       call MPAS_pool_get_config(block % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
       call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
       call MPAS_pool_get_config(block % configs, "config_use_zaerosols",config_use_zaerosols)
       call MPAS_pool_get_config(block % configs, "config_use_vertical_tracers",config_use_vertical_tracers)
       call MPAS_pool_get_config(block % configs, "config_use_atm_dust_file", config_use_atm_dust_file)
       call MPAS_pool_get_config(block % configs, "config_use_iron_solubility_file", config_use_iron_solubility_file)

       call MPAS_pool_get_array(biogeochemistry, "newlyFormedIce", newlyFormedIce)
       call MPAS_pool_get_array(biogeochemistry, "netNitrateUptake", netNitrateUptake)
       call MPAS_pool_get_array(biogeochemistry, "netAmmoniumUptake", netAmmoniumUptake)
       call MPAS_pool_get_array(biogeochemistry, "totalChlorophyll", totalChlorophyll)
       call MPAS_pool_get_array(biogeochemistry, "netSpecificAlgalGrowthRate", netSpecificAlgalGrowthRate)
       call MPAS_pool_get_array(biogeochemistry, "primaryProduction", primaryProduction)
       call MPAS_pool_get_array(biogeochemistry, "netBrineHeight", netBrineHeight)
       call MPAS_pool_get_array(biogeochemistry, "brineBottomChange", brineBottomChange)
       call MPAS_pool_get_array(biogeochemistry, "brineTopChange", brineTopChange)
       call MPAS_pool_get_array(biogeochemistry, "bioPorosity", bioPorosity)
       call MPAS_pool_get_array(biogeochemistry, "bioDiffusivity", bioDiffusivity)
       call MPAS_pool_get_array(biogeochemistry, "bioPermeability", bioPermeability)
       call MPAS_pool_get_array(biogeochemistry, "bioShortwaveFlux", bioShortwaveFlux)
       call MPAS_pool_get_array(biogeochemistry, "bioShortwaveFluxCell", bioShortwaveFluxCell)
       call MPAS_pool_get_array(biogeochemistry, "darcyVelocityBio", darcyVelocityBio)
       call MPAS_pool_get_array(biogeochemistry, "snowIceBioFluxes", snowIceBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "atmosIceBioFluxes", atmosIceBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioConcentrations", oceanBioConcentrations)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioConcentrationsInUse", oceanBioConcentrationsInUse)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioToIceInUse", oceanBioToIceInUse)
       call MPAS_pool_get_array(biogeochemistry, "totalVerticalBiologyIce", totalVerticalBiologyIce)
       call MPAS_pool_get_array(biogeochemistry, "totalVerticalBiologySnow", totalVerticalBiologySnow)
       call MPAS_pool_get_array(biogeochemistry, "atmosBioFluxes", atmosBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "atmosBlackCarbonFlux", atmosBlackCarbonFlux)
       call MPAS_pool_get_array(biogeochemistry, "atmosDustFlux", atmosDustFlux)
       call MPAS_pool_get_array(biogeochemistry, "atmosWetDustFlux", atmosWetDustFlux)
       call MPAS_pool_get_array(biogeochemistry, "atmosDryDustFlux", atmosDryDustFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioFluxes", oceanBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioFluxesCategory", oceanBioFluxesCategory)
       call MPAS_pool_get_array(biogeochemistry, "verticalNitrogenLosses", verticalNitrogenLosses)
       call MPAS_pool_get_array(biogeochemistry, "bioTemperature", bioTemperature)
       call MPAS_pool_get_array(biogeochemistry, "totalSkeletalAlgae", totalSkeletalAlgae)
       call MPAS_pool_get_array(biogeochemistry, "oceanAlgaeConc",oceanAlgaeConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanDOCConc",oceanDOCConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanDICConc",oceanDICConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanDONConc",oceanDONConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanParticulateIronConc",oceanParticulateIronConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanDissolvedIronConc",oceanDissolvedIronConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanNitrateConc",oceanNitrateConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanSilicateConc",oceanSilicateConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanAmmoniumConc",oceanAmmoniumConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSConc",oceanDMSConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSPConc",oceanDMSPConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanHumicsConc",oceanHumicsConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanZAerosolConc",oceanZAerosolConc)

       call MPAS_pool_get_array(diagnostics_biogeochemistry, "bgridPorosityIceCell",bgridPorosityIceCell)
       call MPAS_pool_get_array(diagnostics_biogeochemistry, "bgridSalinityIceCell",bgridSalinityIceCell)
       call MPAS_pool_get_array(diagnostics_biogeochemistry, "bgridTemperatureIceCell",bgridTemperatureIceCell)

       call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTemperature", seaSurfaceTemperature)
       call MPAS_pool_get_array(ocean_coupling, "seaSurfaceSalinity", seaSurfaceSalinity)
       call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)

       call MPAS_pool_get_array(atmos_coupling, "snowfallRate", snowfallRate)

       call MPAS_pool_get_array(atmos_forcing, "IRON_Zolubility_wet", IRON_Zolubility_wet)
       call MPAS_pool_get_array(atmos_forcing, "IRON_Zolubility_dry", IRON_Zolubility_dry)
       call MPAS_pool_get_array(atmos_forcing, "IRON_in_duzt_fraction", IRON_in_duzt_fraction)

       call MPAS_pool_get_array(icestate, "iceAreaCategoryInitial", iceAreaCategoryInitial)
       call MPAS_pool_get_array(icestate, "iceVolumeCategoryInitial", iceVolumeCategoryInitial)
       call MPAS_pool_get_array(icestate, "snowVolumeCategoryInitial", snowVolumeCategoryInitial)
       call MPAS_pool_get_array(icestate, "iceThicknessCategoryInitial", iceThicknessCategoryInitial)
       call MPAS_pool_get_array(icestate, "openWaterArea", openWaterArea)

       call MPAS_pool_get_array(shortwave, "penetratingShortwaveFlux", penetratingShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "shortwaveLayerPenetration", shortwaveLayerPenetration)

       call MPAS_pool_get_array(melt_growth_rates, "basalIceMeltCategory", basalIceMeltCategory)
       call MPAS_pool_get_array(melt_growth_rates, "surfaceIceMeltCategory", surfaceIceMeltCategory)
       call MPAS_pool_get_array(melt_growth_rates, "congelationCategory", congelationCategory)
       call MPAS_pool_get_array(melt_growth_rates, "snowiceFormationCategory", snowiceFormationCategory)
       call MPAS_pool_get_array(melt_growth_rates, "snowMeltCategory", snowMeltCategory)

       call MPAS_pool_get_array(initial,"initialSalinityProfile",initialSalinityProfile)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "skeletalAlgaeConc", skeletalAlgaeConc, 1)
       call MPAS_pool_get_array(tracers, "brineFraction", brineFraction, 1)

       ! newly formed ice
       allocate(newlyFormedIceLogical(nCategories))
       allocate(brineHeightCatInitial(nCategories,nCellsSolve))

       if (checkCarbon) then
          allocate(carbonError(nCellsSolve))
          allocate(totalCarbonCatFinal(nCategories,nCellsSolve))
          allocate(totalCarbonCatInitial(nCategories,nCellsSolve))
          allocate(totalCarbonCatFlux(nCategories,nCellsSolve))
          allocate(brineHeightCatFinal(nCategories,nCellsSolve))
          allocate(totalCarbonFinal(nCellsSolve))
          allocate(totalCarbonInitial(nCellsSolve))
          allocate(totalCarbonFlux(nCellsSolve))
       else
          allocate(carbonError(1))
          allocate(totalCarbonCatFinal(1,1))
          allocate(totalCarbonCatInitial(1,1))
          allocate(totalCarbonCatFlux(1,1))
          allocate(brineHeightCatFinal(1,1))
          allocate(totalCarbonFinal(1))
          allocate(totalCarbonInitial(1))
          allocate(totalCarbonFlux(1))
        endif

       brineHeightCatInitial(:,:) = 0.0_RKIND
       carbonError(:) = 0.0_RKIND
       totalCarbonCatFinal(:,:) = 0.0_RKIND
       totalCarbonCatInitial(:,:) = 0.0_RKIND
       totalCarbonCatFlux(:,:) = 0.0_RKIND
       brineHeightCatFinal(:,:) = 0.0_RKIND
       totalCarbonFinal(:) = 0.0_RKIND
       totalCarbonInitial(:) = 0.0_RKIND
       totalCarbonFlux(:) = 0.0_RKIND

       setGetPhysicsTracers = .true.
       setGetBGCTracers     = (config_use_column_biogeochemistry .or. config_use_zaerosols)

       atmosBioFluxes(:,:) = 0.0_RKIND
       oceanBioConcentrationsInUse(:,:) = 0.0_RKIND
       oceanBioToIceInUse(:,:) = 0.0_RKIND

       !$offomp parallel do default(shared) private(iCategory,iBioTracers,iAlgae, iBioLayers) &
       !$offomp& firstprivate(atmosBioFluxes,atmosBlackCarbonFlux, &
       !$offomp& totalCarbonCatInitial, totalCarbonCatFinal, &
       !$offomp& totalCarbonInitial, totalCarbonFinal, totalCarbonFlux, &
       !$offomp& atmosDustFlux, bioShortwaveFluxCell, newlyFormedIce)
       !
       do iCell = 1, nCellsSolve
          ! newly formed ice
          do iCategory = 1, nCategories
             newlyFormedIceLogical(iCategory) = (newlyFormedIce(iCategory,iCell) == 1)
             brineHeightCatInitial(iCategory,iCell) = brineFraction(1,iCategory,iCell) * &
                iceVolumeCategoryInitial(iCategory,iCell)/(iceAreaCategoryInitial(iCategory,iCell) + seaicePuny)
          enddo ! iCategory

          !update ocean concentrations fields and atmospheric fluxes into allocated array

#ifdef coupled
          call icepack_load_ocean_bio_array(&
                  nit=oceanNitrateConc(iCell), &
                  amm=oceanAmmoniumConc(iCell), &
                  sil=oceanSilicateConc(iCell), &
                  dmsp=oceanDMSPConc(iCell), &
                  dms=oceanDMSConc(iCell), &
                  algalN=oceanAlgaeConc(:,iCell), &
                  doc=oceanDOCConc(:,iCell), &
                  don=oceanDONConc(:,iCell), &
                  dic=oceanDICConc(:,iCell), &
                  fed=oceanDissolvedIronConc(:,iCell), &
                  fep=oceanParticulateIronConc(:,iCell), &
                  zaeros=oceanZAerosolConc(:,iCell), &
                  ocean_bio_all=oceanBioConcentrations(:,iCell), &
                  hum=oceanHumicsConc(iCell))
#else
          do iBioTracers = 1, maxBCType
             atmosBlackCarbonFlux(iBioTracers,iCell) =  1.e-12_RKIND
          enddo
          if (.not. config_use_atm_dust_file) then
             do iBioTracers = 1, maxDustType
                atmosDustFlux(iBioTracers,iCell) =  1.e-13_RKIND
                atmosWetDustFlux(iBioTracers,iCell) =  1.e-13_RKIND * 0.5_RKIND
                atmosDryDustFlux(iBioTracers,iCell) =  1.e-13_RKIND * 0.5_RKIND
             enddo
          endif
#endif
          if (config_use_zaerosols) then
             indexj = ciceTracerObject % index_verticalAerosolsConcLayer(1)
             do iBioTracers = 1, maxBCType
                atmosBioFluxes(indexj -1 + iBioTracers,iCell) = atmosBlackCarbonFlux(iBioTracers,iCell)
             enddo
             do iBioTracers = maxBCType + 1, nzAerosols
                atmosBioFluxes(indexj -1 + iBioTracers, iCell) = atmosDustFlux(iBioTracers-maxBCType,iCell)
             enddo
             indexj = ciceTracerObject % index_dissolvedIronConcLayer(1)
             if (config_use_iron_solubility_file) then
                atmosBioFluxes(indexj,iCell) = 0.0_RKIND
                do iBioTracers = maxBCType + 1, nzAerosols
                   atmosBioFluxes(indexj,iCell) = atmosBioFluxes(indexj,iCell) + &
                      (atmosWetDustFlux(iBioTracers-maxBCType,iCell) * IRON_Zolubility_wet(iCell) + &
                      atmosDryDustFlux(iBioTracers-maxBCType,iCell) * IRON_Zolubility_dry(iCell) ) * &
                      IRON_in_duzt_fraction(iCell) * microgramsPerKilograms / gramsIronPerMolIron
                end do
             end if
          endif

          do iBioTracers = 1, ciceTracerObject % nBioTracers
             iBioData = ciceTracerObject % index_LayerIndexToDataArray(iBioTracers)
             oceanBioConcentrationsInUse(iBioTracers,iCell) =  oceanBioConcentrations(iBioData,iCell)
          enddo ! iBioTracers

          call set_cice_tracer_array_category(block, ciceTracerObject, &
               tracerArrayCategory, iCell, setGetPhysicsTracers, setGetBGCTracers)

          if (checkCarbon) then
             call seaice_total_carbon_content_category(block,&
             totalCarbonCatInitial(:,iCell),iceAreaCategoryInitial(:,:),iceVolumeCategoryInitial(:,:),iCell)
             totalCarbonInitial(iCell) = 0.0_RKIND
             do iCategory = 1, nCategories
                totalCarbonInitial(iCell) = totalCarbonInitial(iCell) + totalCarbonCatInitial(iCategory,iCell)*iceAreaCategoryInitial(iCategory,iCell)
             enddo
            endif

          ! code abort
          abortFlag = .false.
          abortMessage = ""

          call icepack_warnings_clear()
          call icepack_biogeochemistry(&
               dt=config_dt, &
               upNO=netNitrateUptake(iCell), &
               upNH=netAmmoniumUptake(iCell), &
               iDi=bioDiffusivity(:,:,iCell), &
               iki=bioPermeability(:,:,iCell), &
               zfswin=bioShortwaveFlux(:,:,iCell), &
               darcy_V=darcyVelocityBio(:,iCell), &
               grow_net=netSpecificAlgalGrowthRate(iCell), &
               PP_net=primaryProduction(iCell), &
               hbri=netBrineHeight(iCell), &
               dhbr_bot=brineBottomChange(:,iCell), &
               dhbr_top=brineTopChange(:,iCell), &
               Zoo=verticalNitrogenLosses(:,:,iCell), &
               fbio_snoice=snowIceBioFluxes(:,iCell), &
               fbio_atmice=atmosIceBioFluxes(:,iCell), &
               ocean_bio_dh=oceanBioToIceInUse(:,iCell), &
               ocean_bio=oceanBioConcentrationsInUse(:,iCell), &
               first_ice=newlyFormedIceLogical(:), &
               fswpenln=shortwaveLayerPenetration(:,:,iCell), &
               bphi=bioPorosity(:,:,iCell), &
               bTiz=bioTemperature(:,:,iCell), &
               ice_bio_net=totalVerticalBiologyIce(:,iCell), &
               snow_bio_net=totalVerticalBiologySnow(:,iCell), &
               totalChla=totalChlorophyll(iCell), &
               fswthrun=penetratingShortwaveFlux(:,iCell), &
               meltbn=basalIceMeltCategory(:,iCell), &
               melttn=surfaceIceMeltCategory(:,iCell), &
               congeln=congelationCategory(:,iCell), &
               snoicen=snowiceFormationCategory(:,iCell), &
               sst=seaSurfaceTemperature(iCell), &
               sss=seaSurfaceSalinity(iCell), &
               Tf=seaFreezingTemperature(iCell), &
               fsnow=snowfallRate(iCell), &
               meltsn=snowMeltCategory(:,iCell), &
               hin_old=iceThicknessCategoryInitial(:,iCell), &
               flux_bio=oceanBioFluxes(:,iCell), &
               flux_bio_atm=atmosBioFluxes(:,iCell), &
               aicen_init=iceAreaCategoryInitial(:,iCell), &
               vicen_init=iceVolumeCategoryInitial(:,iCell), &
               aicen=iceAreaCategory(1,:,iCell), &
               vicen=iceVolumeCategory(1,:,iCell), &
               vsnon=snowVolumeCategory(1,:,iCell), &
               aice0=openWaterArea(iCell), &
               trcrn=tracerArrayCategory, &
               vsnon_init=snowVolumeCategoryInitial(:,iCell), &
               flux_bion=oceanBioFluxesCategory(:,:,iCell), &
               bioPorosityIceCell=bgridPorosityIceCell(:,iCell), &
               bioSalinityIceCell=bgridSalinityIceCell(:,iCell), &
               bioTemperatureIceCell=bgridTemperatureIceCell(:,iCell))

          abortFlag = icepack_warnings_aborted()

          call get_cice_tracer_array_category(block, ciceTracerObject, &
             tracerArrayCategory, iCell, setGetPhysicsTracers, setGetBGCTracers)

          if (config_use_vertical_tracers) &
             call define_total_biogeochemistry_array_cell(block, ciceTracerObject, totalVerticalBiologyIce(:,iCell), &
             totalVerticalBiologySnow(:,iCell), iCell)

          if (checkCarbon) then
            call seaice_total_carbon_content_category(block,totalCarbonCatFinal(:,iCell),iceAreaCategory(1,:,:),iceVolumeCategory(1,:,:),iCell)
            call seaice_ocean_carbon_flux(block,totalCarbonCatFlux(:,iCell),oceanBioFluxesCategory(:,:,:),iCell)
            totalCarbonFinal(iCell) = 0.0_RKIND
            totalCarbonFlux(iCell) = 0.0_RKIND
            do iCategory = 1, nCategories
               totalCarbonFinal(iCell) = totalCarbonFinal(iCell) + totalCarbonCatFinal(iCategory,iCell)*iceAreaCategory(1,iCategory,iCell)
               totalCarbonFlux(iCell) = totalCarbonFlux(iCell) + totalCarbonCatFlux(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)
            enddo
            carbonError(iCell) = (totalCarbonFinal(iCell) - totalCarbonInitial(iCell))/config_dt +totalCarbonFlux(iCell)
            errorCheck = max(accuracy,accuracy*abs(totalCarbonFlux(iCell)))

            if (abs(carbonError(iCell)) >  errorCheck) then
                  do iCategory = 1,nCategories
                     if (iceAreaCategory(1,iCategory,iCell) > seaicePuny) then
                     brineHeightCatFinal(iCategory,iCell) = brineFraction(1,iCategory,iCell) * &
                     iceVolumeCategory(1,iCategory,iCell)/(iceAreaCategory(1,iCategory,iCell) + seaicePuny)
                     call mpas_log_write("column_biogeochemistry, carbon conservation error", messageType=MPAS_LOG_ERR)
                     call mpas_log_write("iCell: $i", messageType=MPAS_LOG_ERR, intArgs=(/indexToCellID(iCell)/))
                     call mpas_log_write("iCategory: $i", messageType=MPAS_LOG_ERR, intArgs=(/iCategory/))
                     call mpas_log_write("carbonError: $r", messageType=MPAS_LOG_ERR, realArgs=(/carbonError(iCell)/))
                     call mpas_log_write("carbonError*iceAreaCategory: $r", messageType=MPAS_LOG_ERR, realArgs=(/carbonError(iCell)*iceAreaCategory(1,iCategory,iCell)/))
                     call mpas_log_write("iceAreaCategory: $r", messageType=MPAS_LOG_ERR, realArgs=(/iceAreaCategory(1,iCategory,iCell)/))
                     call mpas_log_write("iceAreaCategoryInitial: $r", messageType=MPAS_LOG_ERR, realArgs=(/iceAreaCategoryInitial(iCategory,iCell)/))
                     call mpas_log_write("totalCarbonCatInitial(iCategory): $r", messageType=MPAS_LOG_ERR, realArgs=(/totalCarbonCatInitial(iCategory,iCell)/))
                     call mpas_log_write("totalCarbonCatFinal(iCategory,iCell): $r", messageType=MPAS_LOG_ERR, realArgs=(/totalCarbonCatFinal(iCategory,iCell)/))
                     call mpas_log_write("totalCarbonCatFlux(iCategory,iCell): $r", messageType=MPAS_LOG_ERR, realArgs=(/totalCarbonCatFlux(iCategory,iCell)/))
                     call mpas_log_write("brineHeightCatInitial(iCategory,iCell): $r", messageType=MPAS_LOG_ERR, realArgs=(/brineHeightCatInitial(iCategory,iCell)/))
                     call mpas_log_write("brineHeightCatFinal(iCategory,iCell): $r", messageType=MPAS_LOG_ERR, realArgs=(/brineHeightCatFinal(iCategory,iCell)/))
                     call mpas_log_write("iceAreaCategory(1,iCategory,iCell): $r", messageType=MPAS_LOG_ERR, realArgs=(/iceAreaCategory(1,iCategory,iCell)/))
                  end if
               enddo            !categories
            endif               ! carbonError
         endif                  ! checkCarbon

          totalSkeletalAlgae(iCell) = 0.0_RKIND
          bioShortwaveFluxCell(:,iCell) = 0.0_RKIND

          do iCategory = 1, nCategories
             if (config_use_skeletal_biochemistry .and. iceAreaCategory(1,iCategory,iCell) > seaicePuny) then
                do iAlgae = 1, nAlgae
                   totalSkeletalAlgae(iCell) = totalSkeletalAlgae(iCell) + &
                        skeletalAlgaeConc(iAlgae,iCategory,iCell) * &
                        iceAreaCategory(1,iCategory,iCell)
                enddo
             endif
             if (config_use_vertical_tracers) then
                do iBioLayers = 1, nBioLayersP1
                   bioShortwaveFluxCell(iBioLayers,iCell) = bioShortwaveFluxCell(iBioLayers,iCell) + &
                        bioShortwaveFlux(iBioLayers,iCategory,iCell) * &
                        iceAreaCategory(1,iCategory,iCell)
                enddo
             endif
             newlyFormedIce(iCategory,iCell) = 0
             if (newlyFormedIceLogical(iCategory)) newlyFormedIce(iCategory,iCell) = 1
          enddo ! iCategory

          ! code abort
          if (abortFlag) then
             call mpas_log_write("column_biogeochemistry: "//trim(abortMessage) , messageType=MPAS_LOG_ERR)
             call mpas_log_write("iCell: $i", messageType=MPAS_LOG_ERR, intArgs=(/indexToCellID(iCell)/))
          endif

          call seaice_icepack_write_warnings(abortFlag)

       enddo ! iCell

       call seaice_critical_error_write_block(domain, block, abortFlag)
       call seaice_check_critical_error(domain, abortFlag)

       deallocate(totalCarbonCatFinal)
       deallocate(totalCarbonCatInitial)
       deallocate(totalCarbonCatFlux)
       deallocate(brineHeightCatFinal)
       deallocate(totalCarbonFinal)
       deallocate(totalCarbonInitial)
       deallocate(totalCarbonFlux)
       deallocate(carbonError)

       deallocate(brineHeightCatInitial)
       deallocate(newlyFormedIceLogical)

       block => block % next
    end do

  end subroutine column_biogeochemistry

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

    type(MPAS_Time_type) :: &
         currentTime

    integer :: &
         dayOfYearInt, &
         ierr

    currentTime = MPAS_get_clock_time(clock, MPAS_NOW, ierr=ierr)

    call MPAS_get_time(currentTime, DoY=dayOfYearInt, ierr=ierr)

    dayOfYear = real(dayOfYearInt, RKIND)

  end subroutine get_day_of_year

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_seconds_into_day
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 4th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_seconds_into_day(clock, secondsIntoDay)

    type(MPAS_clock_type), intent(in) :: &
         clock

    integer, intent(out) :: &
         secondsIntoDay

    type(MPAS_Time_type) :: &
         currentTime

    integer :: &
         ierr, &
         hours, &
         minutes, &
         seconds

    currentTime = MPAS_get_clock_time(clock, MPAS_NOW, ierr=ierr)

    call MPAS_get_time(currentTime, H=hours, M=minutes, S=seconds, ierr=ierr)

    secondsIntoDay = hours * 3600 + minutes * 60 + seconds

  end subroutine get_seconds_into_day

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_days_in_year
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 4th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_days_in_year(domain, clock, daysInYear)

    type(domain_type), intent(in) :: domain

    type(MPAS_clock_type), intent(in) :: &
         clock

    integer, intent(out) :: &
         daysInYear

    type(MPAS_Time_type) :: &
         currentTime

    character(len=strKIND), pointer :: &
         config_calendar_type

    integer :: &
         ierr, &
         year

    currentTime = MPAS_get_clock_time(clock, MPAS_NOW, ierr=ierr)

    call MPAS_get_time(currentTime, YYYY=year, ierr=ierr)

    call MPAS_pool_get_config(domain % configs, "config_calendar_type", config_calendar_type)

    select case (trim(config_calendar_type))
    case ("gregorian")
       if (isLeapYear(Year)) then
          daysInYear = sum(daysInMonthLeap)
       else
          daysInYear = sum(daysInMonth)
       endif
    case ("noleap")
       daysInYear = sum(daysInMonth)
    end select

  end subroutine get_days_in_year

!-----------------------------------------------------------------------
! Other routines
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_column_update_state
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_column_update_state(domain, stateUpdateType, dt, iceAgeTimeOffset)

    type(domain_type), intent(inout) :: domain

    character(len=*), intent(in) :: &
         stateUpdateType

    real(kind=RKIND), intent(in) :: &
         dt, &
         iceAgeTimeOffset

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         tracers_aggregate, &
         diagnostics

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         iceVolumeCell, &
         iceAgeCell, &
         iceAreaTendency, &
         iceVolumeTendency, &
         iceAgeTendency

    integer, pointer :: &
         nCellsSolve

    integer :: &
         iCell

    logical, pointer :: &
         config_use_ice_age

    ! aggregate state variables
    call mpas_timer_start("Column aggregate")
    call seaice_icepack_aggregate(domain)
    call mpas_timer_stop("Column aggregate")

    ! get configs
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_ice_age", config_use_ice_age)

    ! compute thermodynamic area and volume tendencies
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
       call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnostics)

       call MPAS_pool_get_dimension(tracers_aggregate, "nCellsSolve", nCellsSolve)

       call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(tracers_aggregate, "iceVolumeCell", iceVolumeCell)
       call MPAS_pool_get_array(tracers_aggregate, "iceAgeCell", iceAgeCell)

       if (trim(stateUpdateType) == "transport") then

          call MPAS_pool_get_array(diagnostics, "iceAreaTendencyTransport", iceAreaTendency)
          call MPAS_pool_get_array(diagnostics, "iceVolumeTendencyTransport", iceVolumeTendency)
          call MPAS_pool_get_array(diagnostics, "iceAgeTendencyTransport", iceAgeTendency)

       else if (trim(stateUpdateType) == "thermodynamics") then

          call MPAS_pool_get_array(diagnostics, "iceAreaTendencyThermodynamics", iceAreaTendency)
          call MPAS_pool_get_array(diagnostics, "iceVolumeTendencyThermodynamics", iceVolumeTendency)
          call MPAS_pool_get_array(diagnostics, "iceAgeTendencyThermodynamics", iceAgeTendency)

       else

          call mpas_log_write("seaice_column_update_state: Unknown update type: "//trim(stateUpdateType), messageType=MPAS_LOG_CRIT)

       endif

       do iCell = 1, nCellsSolve

          iceAreaTendency(iCell)   = (iceAreaCell(iCell)   - iceAreaTendency(iCell))   / dt
          iceVolumeTendency(iCell) = (iceVolumeCell(iCell) - iceVolumeTendency(iCell)) / dt

          if (config_use_ice_age) then
             if (iceAgeTimeOffset > 0.0_RKIND) then

                if (iceAgeCell(iCell) > 0.0_RKIND) &
                     iceAgeTendency(iCell) = &
                          (iceAgeCell(iCell) - iceAgeTendency(iCell) - iceAgeTimeOffset) / dt

             else

                iceAgeTendency(iCell) = &
                     (iceAgeCell(iCell) - iceAgeTendency(iCell)) / dt

             endif
          endif

       enddo ! iCell

       block => block % next
    enddo

  end subroutine seaice_column_update_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack_aggregate
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_icepack_aggregate(domain)

    use icepack_intfc, only: icepack_aggregate

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate, &
         icestate, &
         ocean_coupling

    logical, pointer :: &
         config_use_column_biogeochemistry, &
         config_use_zaerosols

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         openWaterArea, &
         seaFreezingTemperature

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    integer :: &
         iCell

    integer, pointer :: &
         nCellsSolve

    integer, dimension(:), pointer :: &
         indexToCellID

    logical :: &
         setGetPhysicsTracers, &
         setGetBGCTracers

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

       call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
       call MPAS_pool_get_config(block % configs, "config_use_zaerosols", config_use_zaerosols)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_array(mesh, "indexToCellID", indexToCellID)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

       call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(tracers_aggregate, "iceVolumeCell", iceVolumeCell)
       call MPAS_pool_get_array(tracers_aggregate, "snowVolumeCell", snowVolumeCell)

       call MPAS_pool_get_array(icestate, "openWaterArea", openWaterArea)

       call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)

       setGetPhysicsTracers = .true.
       setGetBGCTracers     = (config_use_column_biogeochemistry .or. config_use_zaerosols)

       do iCell = 1, nCellsSolve

          ! set the category tracer array
          call set_cice_tracer_array_category(block, ciceTracerObject, &
               tracerArrayCategory, iCell, setGetPhysicsTracers, setGetBGCTracers)

          call icepack_aggregate(                    &
               iceAreaCategory(1,:,iCell),           &
               tracerArrayCategory,                  & ! trcrn
               iceVolumeCategory(1,:,iCell),         &
               snowVolumeCategory(1,:,iCell),        &
               iceAreaCell(iCell),                   &
               tracerArrayCell,                      & ! trcr
               iceVolumeCell(iCell),                 &
               snowVolumeCell(iCell),                &
               openWaterArea(iCell),                 &
               ciceTracerObject % parentIndex,       & ! trcr_depend
               ciceTracerObject % firstAncestorMask, & ! trcr_base
               ciceTracerObject % ancestorNumber,    & ! n_trcr_strata
               ciceTracerObject % ancestorIndices,   & ! nt_strata
               seaFreezingTemperature(iCell))          ! Tf

          ! set the cell tracer array
          call get_cice_tracer_array_cell(block, ciceTracerObject, &
               tracerArrayCell, iCell, setGetPhysicsTracers, setGetBGCTracers)

       enddo ! iCell

       block => block % next
    end do

    call seaice_icepack_write_warnings(icepack_warnings_aborted())

  end subroutine seaice_icepack_aggregate

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack_aggregate_simple
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_icepack_aggregate_simple(domain)

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate, &
         icestate

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         openWaterArea

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    integer :: &
         iCell, &
         iCategory

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

       call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(tracers_aggregate, "iceVolumeCell", iceVolumeCell)
       call MPAS_pool_get_array(tracers_aggregate, "snowVolumeCell", snowVolumeCell)

       call MPAS_pool_get_array(icestate, "openWaterArea", openWaterArea)

       do iCell = 1, nCellsSolve

          iceAreaCell(iCell)   = 0.0_RKIND
          iceVolumeCell(iCell) = 0.0_RKIND

          do iCategory = 1, nCategories

             iceAreaCell(iCell)   = iceAreaCell(iCell)   + iceAreaCategory(1,iCategory,iCell)
             iceVolumeCell(iCell) = iceVolumeCell(iCell) + iceVolumeCategory(1,iCategory,iCell)

          enddo ! iCategory

          openWaterArea(iCell) = 1.0_RKIND - iceAreaCell(iCell)

       enddo ! iCell

       block => block % next
    end do

  end subroutine seaice_icepack_aggregate_simple

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack_coupling_prep
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_icepack_coupling_prep(domain)

    use seaice_constants, only: &
         seaicePuny, &
         seaiceDensityFreshwater

    type(domain_type) :: domain

    type(block_type), pointer :: block

    logical, pointer :: &
         config_use_ocean_mixed_layer, &
         config_include_pond_freshwater_feedback, &
         config_use_column_biogeochemistry, &
         config_use_zaerosols

    type(MPAS_pool_type), pointer :: &
         oceanCoupling, &
         diagnostics, &
         shortwave, &
         atmosCoupling, &
         tracers, &
         ponds, &
         oceanFluxes, &
         biogeochemistry, &
         mesh

    real(kind=RKIND), dimension(:), pointer :: &
         freezingMeltingPotential, &
         freezingMeltingPotentialInitial, &
         albedoVisibleDirectCell, &
         albedoVisibleDiffuseCell, &
         albedoIRDirectCell, &
         albedoIRDiffuseCell, &
         albedoVisibleDirectArea, &
         albedoVisibleDiffuseArea, &
         albedoIRDirectArea, &
         albedoIRDiffuseArea, &
         bareIceAlbedoCell, &
         snowAlbedoCell, &
         pondAlbedoCell, &
         solarZenithAngleCosine, &
         effectivePondAreaCell, &
         shortwaveScalingFactor, &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown, &
         pondFreshWaterFlux, &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         oceanShortwaveFlux, &
         oceanFreshWaterFluxArea, &
         oceanSaltFluxArea, &
         oceanHeatFluxArea, &
         oceanShortwaveFluxArea, &
         oceanNitrateFlux, &
         oceanSilicateFlux, &
         oceanAmmoniumFlux, &
         oceanDMSFlux, &
         oceanDMSPpFlux, &
         oceanDMSPdFlux, &
         oceanHumicsFlux, &
         oceanDustIronFlux, &
         oceanBlackCarbonFlux, &
         totalOceanCarbonFlux, &
         oceanNitrateFluxArea, &
         oceanSilicateFluxArea, &
         oceanAmmoniumFluxArea, &
         oceanDMSFluxArea, &
         oceanDMSPpFluxArea, &
         oceanDMSPdFluxArea, &
         oceanHumicsFluxArea, &
         oceanDustIronFluxArea, &
         oceanBlackCarbonFluxArea

    real(kind=RKIND), dimension(:,:), pointer :: &
         albedoVisibleDirectCategory, &
         albedoVisibleDiffuseCategory, &
         albedoIRDirectCategory, &
         albedoIRDiffuseCategory, &
         bareIceAlbedoCategory, &
         snowAlbedoCategory, &
         pondAlbedoCategory, &
         effectivePondAreaCategory, &
         oceanBioFluxes, &
         oceanAlgaeFlux, &
         oceanDOCFlux, &
         oceanDICFlux, &
         oceanDONFlux, &
         oceanParticulateIronFlux, &
         oceanDissolvedIronFlux, &
         oceanAlgaeFluxArea, &
         oceanDOCFluxArea, &
         oceanDICFluxArea, &
         oceanDONFluxArea, &
         oceanParticulateIronFluxArea, &
         oceanDissolvedIronFluxArea

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory

    real(kind=RKIND), pointer :: &
         config_dt

    real(kind=RKIND), pointer :: &
         config_ratio_C_to_N_diatoms, &
         config_ratio_C_to_N_small_plankton, &
         config_ratio_C_to_N_phaeocystis, &
         config_ratio_C_to_N_proteins

    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nZBGCTracers, &
         maxAlgaeType, &
         maxDOCType, &
         maxDICType, &
         maxDONType, &
         maxIronType, &
         maxBCType, &
         maxDustType, &
         maxAerosolType

    integer :: &
         iCell, &
         iCategory, &
         iBioTracers, &
         iBioData

    real(kind=RKIND), dimension(:), allocatable :: &
         ratio_C_to_N

    real(kind=RKIND), dimension(:), allocatable :: &
         oceanBioFluxesAll

    call MPAS_pool_get_config(domain % configs, "config_use_ocean_mixed_layer", config_use_ocean_mixed_layer)
    call MPAS_pool_get_config(domain % configs, "config_dt", config_dt)
    call MPAS_pool_get_config(domain % configs, "config_include_pond_freshwater_feedback", config_include_pond_freshwater_feedback)
    call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_zaerosols", config_use_zaerosols)
    call MPAS_pool_get_config(domain % configs, "config_ratio_C_to_N_diatoms", config_ratio_C_to_N_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_ratio_C_to_N_small_plankton", config_ratio_C_to_N_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_ratio_C_to_N_phaeocystis", config_ratio_C_to_N_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_ratio_C_to_N_proteins", config_ratio_C_to_N_proteins)

    if (config_use_ocean_mixed_layer) &
         call seaice_icepack_ocean_mixed_layer(domain)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCoupling)
       call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnostics)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
       call MPAS_pool_get_subpool(block % structs, "ponds", ponds)
       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxes)
       call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)
       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

       call MPAS_pool_get_dimension(tracers, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(tracers, "nCategories", nCategories)

       call MPAS_pool_get_array(oceanCoupling, "freezingMeltingPotential", freezingMeltingPotential)
       call MPAS_pool_get_array(diagnostics, "freezingMeltingPotentialInitial", freezingMeltingPotentialInitial)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)

       call MPAS_pool_get_array(shortwave, "albedoVisibleDirectCell", albedoVisibleDirectCell)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDiffuseCell", albedoVisibleDiffuseCell)
       call MPAS_pool_get_array(shortwave, "albedoIRDirectCell", albedoIRDirectCell)
       call MPAS_pool_get_array(shortwave, "albedoIRDiffuseCell", albedoIRDiffuseCell)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDirectCategory", albedoVisibleDirectCategory)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDiffuseCategory", albedoVisibleDiffuseCategory)
       call MPAS_pool_get_array(shortwave, "albedoIRDirectCategory", albedoIRDirectCategory)
       call MPAS_pool_get_array(shortwave, "albedoIRDiffuseCategory", albedoIRDiffuseCategory)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDirectArea", albedoVisibleDirectArea)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDiffuseArea", albedoVisibleDiffuseArea)
       call MPAS_pool_get_array(shortwave, "albedoIRDirectArea", albedoIRDirectArea)
       call MPAS_pool_get_array(shortwave, "albedoIRDiffuseArea", albedoIRDiffuseArea)

       call MPAS_pool_get_array(shortwave, "solarZenithAngleCosine", solarZenithAngleCosine)
       call MPAS_pool_get_array(shortwave, "bareIceAlbedoCell", bareIceAlbedoCell)
       call MPAS_pool_get_array(shortwave, "snowAlbedoCell", snowAlbedoCell)
       call MPAS_pool_get_array(shortwave, "pondAlbedoCell", pondAlbedoCell)
       call MPAS_pool_get_array(shortwave, "bareIceAlbedoCategory", bareIceAlbedoCategory)
       call MPAS_pool_get_array(shortwave, "snowAlbedoCategory", snowAlbedoCategory)
       call MPAS_pool_get_array(shortwave, "pondAlbedoCategory", pondAlbedoCategory)

       call MPAS_pool_get_array(shortwave, "effectivePondAreaCell", effectivePondAreaCell)
       call MPAS_pool_get_array(shortwave, "effectivePondAreaCategory", effectivePondAreaCategory)

       call MPAS_pool_get_array(shortwave, "shortwaveScalingFactor", shortwaveScalingFactor)

       call MPAS_pool_get_array(atmosCoupling, "shortwaveVisibleDirectDown", shortwaveVisibleDirectDown)
       call MPAS_pool_get_array(atmosCoupling, "shortwaveVisibleDiffuseDown", shortwaveVisibleDiffuseDown)
       call MPAS_pool_get_array(atmosCoupling, "shortwaveIRDirectDown", shortwaveIRDirectDown)
       call MPAS_pool_get_array(atmosCoupling, "shortwaveIRDiffuseDown", shortwaveIRDiffuseDown)

       call MPAS_pool_get_array(ponds, "pondFreshWaterFlux", pondFreshWaterFlux)

       call MPAS_pool_get_array(oceanFluxes, "oceanFreshWaterFlux", oceanFreshWaterFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanSaltFlux", oceanSaltFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanHeatFlux", oceanHeatFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanShortwaveFlux", oceanShortwaveFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanFreshWaterFluxArea", oceanFreshWaterFluxArea)
       call MPAS_pool_get_array(oceanFluxes, "oceanSaltFluxArea", oceanSaltFluxArea)
       call MPAS_pool_get_array(oceanFluxes, "oceanHeatFluxArea", oceanHeatFluxArea)
       call MPAS_pool_get_array(oceanFluxes, "oceanShortwaveFluxArea", oceanShortwaveFluxArea)

       call MPAS_pool_get_array(biogeochemistry, "oceanNitrateFlux", oceanNitrateFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanSilicateFlux", oceanSilicateFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanAmmoniumFlux", oceanAmmoniumFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSFlux", oceanDMSFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSPpFlux", oceanDMSPpFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSPdFlux", oceanDMSPdFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanHumicsFlux", oceanHumicsFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDustIronFlux", oceanDustIronFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanBlackCarbonFlux", oceanBlackCarbonFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioFluxes", oceanBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "oceanAlgaeFlux", oceanAlgaeFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDOCFlux", oceanDOCFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDICFlux", oceanDICFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDONFlux", oceanDONFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanParticulateIronFlux", oceanParticulateIronFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDissolvedIronFlux", oceanDissolvedIronFlux)
       call MPAS_pool_get_array(biogeochemistry, "totalOceanCarbonFlux", totalOceanCarbonFlux)

       call MPAS_pool_get_array(biogeochemistry, "oceanNitrateFluxArea", oceanNitrateFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanSilicateFluxArea", oceanSilicateFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanAmmoniumFluxArea", oceanAmmoniumFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSFluxArea", oceanDMSFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSPpFluxArea", oceanDMSPpFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSPdFluxArea", oceanDMSPdFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanHumicsFluxArea", oceanHumicsFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanDustIronFluxArea", oceanDustIronFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanBlackCarbonFluxArea", oceanBlackCarbonFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanAlgaeFluxArea", oceanAlgaeFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanDOCFluxArea", oceanDOCFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanDICFluxArea", oceanDICFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanDONFluxArea", oceanDONFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanParticulateIronFluxArea", oceanParticulateIronFluxArea)
       call MPAS_pool_get_array(biogeochemistry, "oceanDissolvedIronFluxArea", oceanDissolvedIronFluxArea)

       call MPAS_pool_get_dimension(mesh, "nZBGCTracers", nZBGCTracers)
       call MPAS_pool_get_dimension(mesh, "maxAlgaeType", maxAlgaeType)
       call MPAS_pool_get_dimension(mesh, "maxDOCType", maxDOCType)
       call MPAS_pool_get_dimension(mesh, "maxDICType", maxDICType)
       call MPAS_pool_get_dimension(mesh, "maxDONType", maxDONType)
       call MPAS_pool_get_dimension(mesh, "maxAerosolType", maxAerosolType)
       call MPAS_pool_get_dimension(mesh, "maxIronType", maxIronType)
       call MPAS_pool_get_dimension(mesh, "maxBCType", maxBCType)
       call MPAS_pool_get_dimension(mesh, "maxDustType", maxDustType)

       allocate(oceanBioFluxesAll(nZBGCTracers))

       allocate(ratio_C_to_N(3))

       ratio_C_to_N(1) = config_ratio_C_to_N_diatoms
       ratio_C_to_N(2) = config_ratio_C_to_N_small_plankton
       ratio_C_to_N(3) = config_ratio_C_to_N_phaeocystis

       do iCell = 1, nCellsSolve

          !-------------------------------------------------------------------
          ! store initial freezing melting potential
          !-------------------------------------------------------------------

          freezingMeltingPotentialInitial(iCell) = freezingMeltingPotential(iCell)

          !-------------------------------------------------------------------
          ! aggregate albedos
          !-------------------------------------------------------------------

          albedoVisibleDirectCell(iCell)  = 0.0_RKIND
          albedoVisibleDiffuseCell(iCell) = 0.0_RKIND
          albedoIRDirectCell(iCell)       = 0.0_RKIND
          albedoIRDiffuseCell(iCell)      = 0.0_RKIND

          bareIceAlbedoCell(iCell) = 0.0_RKIND
          snowAlbedoCell(iCell)    = 0.0_RKIND
          pondAlbedoCell(iCell)    = 0.0_RKIND

          effectivePondAreaCell(iCell) = 0.0_RKIND

          do iCategory = 1, nCategories

             albedoVisibleDirectCell(iCell) = albedoVisibleDirectCell(iCell) + &
                  albedoVisibleDirectCategory(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)
             albedoVisibleDiffuseCell(iCell) = albedoVisibleDiffuseCell(iCell) + &
                  albedoVisibleDiffuseCategory(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)
             albedoIRDirectCell(iCell) = albedoIRDirectCell(iCell) + &
                  albedoIRDirectCategory(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)
             albedoIRDiffuseCell(iCell) = albedoIRDiffuseCell(iCell) + &
                  albedoIRDiffuseCategory(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)

             ! sun above horizon
             if (solarZenithAngleCosine(iCell) > seaicePuny) then

                bareIceAlbedoCell(iCell) = bareIceAlbedoCell(iCell) + &
                     bareIceAlbedoCategory(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)
                snowAlbedoCell(iCell) = snowAlbedoCell(iCell) + &
                     snowAlbedoCategory(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)
                pondAlbedoCell(iCell) = pondAlbedoCell(iCell) + &
                     pondAlbedoCategory(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)

             endif

             effectivePondAreaCell(iCell) = effectivePondAreaCell(iCell) + &
                  effectivePondAreaCategory(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)

          enddo ! iCategory

          !-------------------------------------------------------------------
          ! reduce oceanFreshWaterFlux by pondFreshWaterFlux for coupling
          !-------------------------------------------------------------------

          if (config_include_pond_freshwater_feedback) then
             pondFreshWaterFlux(iCell)  = pondFreshWaterFlux(iCell) * seaiceDensityFreshwater / config_dt
             oceanFreshWaterFlux(iCell) = oceanFreshWaterFlux(iCell) - pondFreshWaterFlux(iCell)
          endif

          !-------------------------------------------------------------------
          ! Store grid box mean albedos and fluxes before scaling by aice
          !-------------------------------------------------------------------

          albedoVisibleDirectArea(iCell)  = albedoVisibleDirectCell(iCell)
          albedoVisibleDiffuseArea(iCell) = albedoVisibleDiffuseCell(iCell)
          albedoIRDirectArea(iCell)       = albedoIRDirectCell(iCell)
          albedoIRDiffuseArea(iCell)      = albedoIRDiffuseCell(iCell)
          oceanFreshWaterFluxArea(iCell)  = oceanFreshWaterFlux(iCell)
          oceanSaltFluxArea(iCell)        = oceanSaltFlux(iCell)
          oceanHeatFluxArea(iCell)        = oceanHeatFlux(iCell)
          oceanShortwaveFluxArea(iCell)   = oceanShortwaveFlux(iCell)

          !-----------------------------------------------------------------
          ! Save net shortwave for scaling factor in shortwaveScalingFactor
          !-----------------------------------------------------------------

          shortwaveScalingFactor(iCell) = &
               shortwaveVisibleDirectDown(iCell)  * (1.0_RKIND - albedoVisibleDirectArea(iCell)) + &
               shortwaveVisibleDiffuseDown(iCell) * (1.0_RKIND - albedoVisibleDiffuseArea(iCell)) + &
               shortwaveIRDirectDown(iCell)       * (1.0_RKIND - albedoIRDirectArea(iCell)) + &
               shortwaveIRDiffuseDown(iCell)      * (1.0_RKIND - albedoIRDiffuseArea(iCell))

          !-----------------------------------------------------------------
          ! Define ocean biogeochemical flux variables
          !-----------------------------------------------------------------

           oceanBioFluxesAll(:)              = 0.0_RKIND

           if (config_use_column_biogeochemistry .or. config_use_zaerosols) then

              totalOceanCarbonFlux(iCell)       = 0.0_RKIND
              oceanAlgaeFlux(:,iCell)           = 0.0_RKIND
              oceanDOCFlux(:,iCell)             = 0.0_RKIND
              oceanDICFlux(:,iCell)             = 0.0_RKIND
              oceanDONFlux(:,iCell)             = 0.0_RKIND
              oceanParticulateIronFlux(:,iCell) = 0.0_RKIND
              oceanDissolvedIronFlux(:,iCell)   = 0.0_RKIND
              oceanNitrateFlux(iCell)           = 0.0_RKIND
              oceanSilicateFlux(iCell)          = 0.0_RKIND
              oceanAmmoniumFlux(iCell)          = 0.0_RKIND
              oceanDMSPpFlux(iCell)             = 0.0_RKIND
              oceanDMSPdFlux(iCell)             = 0.0_RKIND
              oceanDMSFlux(iCell)               = 0.0_RKIND
              oceanDustIronFlux(iCell)          = 0.0_RKIND
              oceanBlackCarbonFlux(iCell)       = 0.0_RKIND
              oceanHumicsFlux(iCell)            = 0.0_RKIND

              do iBioTracers = 1, ciceTracerObject % nBioTracers
                 iBioData = ciceTracerObject % index_LayerIndexToDataArray(iBioTracers)
                 oceanBioFluxesAll(iBioData) = oceanBioFluxes(iBioTracers,iCell)
              enddo
              iBioData = 0

              ! Algae
              do iBioTracers = 1, maxAlgaeType
                  iBioData = iBioData+1
                  oceanAlgaeFlux(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
                  oceanAlgaeFluxArea(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
                  totalOceanCarbonFlux(iCell) = totalOceanCarbonFlux(iCell) + &
                      oceanAlgaeFlux(iBioTracers,iCell) * ratio_C_to_N(iBioTracers)
              enddo

              ! Nitrate
              iBioData = iBioData+1
              oceanNitrateFlux(iCell) = oceanBioFluxesAll(iBioData)
              oceanNitrateFluxArea(iCell) = oceanBioFluxesAll(iBioData)

              ! Polysaccharids and Lipids
              do iBioTracers = 1, maxDOCType
                  iBioData = iBioData+1
                  oceanDOCFlux(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
                  oceanDOCFluxArea(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
                  totalOceanCarbonFlux(iCell) = totalOceanCarbonFlux(iCell) + &
                      oceanDOCFlux(iBioTracers,iCell)
              enddo

              ! DIC
              do iBioTracers = 1, maxDICType
                  iBioData = iBioData+1
                  oceanDICFlux(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
                  oceanDICFluxArea(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
                  totalOceanCarbonFlux(iCell) = totalOceanCarbonFlux(iCell) + &
                      oceanDICFlux(iBioTracers,iCell)
              enddo

              ! Chlorophyll (not saved)
              iBioData = iBioData+maxAlgaeType

              ! Ammonium
              iBioData = iBioData+1
              oceanAmmoniumFlux(iCell) = oceanBioFluxesAll(iBioData)
              oceanAmmoniumFluxArea(iCell) = oceanBioFluxesAll(iBioData)

              ! Silicate
              iBioData = iBioData+1
              oceanSilicateFlux(iCell) = oceanBioFluxesAll(iBioData)
              oceanSilicateFluxArea(iCell) = oceanBioFluxesAll(iBioData)

              ! DMSPp
              iBioData = iBioData+1
              oceanDMSPpFlux(iCell) = oceanBioFluxesAll(iBioData)
              oceanDMSPpFluxArea(iCell) = oceanBioFluxesAll(iBioData)

              ! DMSPd
              iBioData = iBioData+1
              oceanDMSPdFlux(iCell) = oceanBioFluxesAll(iBioData)
              oceanDMSPdFluxArea(iCell) = oceanBioFluxesAll(iBioData)

              ! DMS
              iBioData = iBioData+1
              oceanDMSFlux(iCell) = oceanBioFluxesAll(iBioData)
              oceanDMSFluxArea(iCell) = oceanBioFluxesAll(iBioData)

              ! PON
              iBioData = iBioData+1

              ! DON (Proteins)
              do iBioTracers = 1, maxDONType
                  iBioData = iBioData+1
                  oceanDONFlux(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
                  oceanDONFluxArea(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
                  totalOceanCarbonFlux(iCell) = totalOceanCarbonFlux(iCell) + &
                      oceanDONFlux(iBioTracers,iCell) * config_ratio_C_to_N_proteins
              enddo

              ! Dissolved Iron
              do iBioTracers = 1, maxIronType
                  iBioData = iBioData+1
                  oceanDissolvedIronFlux(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
                  oceanDissolvedIronFluxArea(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
              enddo

              ! Particulate Iron
              do iBioTracers = 1, maxIronType
                  iBioData = iBioData+1
                  oceanParticulateIronFlux(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
                  oceanParticulateIronFluxArea(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
              enddo

              ! Black Carbon (combined; saved for conservation)
              do iBioTracers = 1, maxBCType
                  iBioData = iBioData + 1
                  oceanBlackCarbonFlux(iCell) = oceanBlackCarbonFlux(iCell) +  oceanBioFluxesAll(iBioData)
               enddo
               oceanBlackCarbonFluxArea(iCell) = oceanBlackCarbonFlux(iCell)

              ! Dust (combined)
              do iBioTracers = 1, maxDustType
                  iBioData = iBioData+1
                  oceanDustIronFlux(iCell) = oceanDustIronFlux(iCell) +  oceanBioFluxesAll(iBioData)
               enddo
               oceanDustIronFluxArea(iCell) = oceanDustIronFlux(iCell)

              ! Humics
              iBioData = iBioData+1
              oceanHumicsFlux(iCell) = oceanBioFluxesAll(iBioData)
              oceanHumicsFluxArea(iCell) = oceanBioFluxesAll(iBioData)
              totalOceanCarbonFlux(iCell) = totalOceanCarbonFlux(iCell) + &
                 oceanHumicsFlux(iCell)

          endif ! config_use_column_biogeochemistry .or. config_use_zaerosols

       enddo ! iCell

       deallocate(oceanBioFluxesAll)
       deallocate(ratio_C_to_N)

       block => block % next
    enddo

    call seaice_column_scale_fluxes(domain)

  end subroutine seaice_icepack_coupling_prep

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_column_scale_fluxes
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_column_scale_fluxes(domain)

    use seaice_constants, only: &
         seaiceStefanBoltzmann, &
         seaiceFreshWaterFreezingPoint

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         tracersAggregate, &
         velocitySolver, &
         atmosFluxes, &
         shortwave, &
         atmosCoupling, &
         oceanCoupling, &
         oceanFluxes, &
         biogeochemistry, &
         mesh

    logical, pointer :: &
         config_use_column_biogeochemistry, &
         config_use_zaerosols

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         airStressCellU, &
         airStressCellV, &
         sensibleHeatFlux, &
         latentHeatFlux, &
         absorbedShortwaveFlux, &
         longwaveUp, &
         evaporativeWaterFlux, &
         atmosReferenceHumidity2m, &
         atmosReferenceTemperature2m, &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         oceanShortwaveFlux, &
         albedoVisibleDirectCell, &
         albedoIRDirectCell, &
         albedoVisibleDiffuseCell, &
         albedoIRDiffuseCell, &
         airTemperature, &
         airSpecificHumidity, &
         seaFreezingTemperature, &
         oceanNitrateFlux, &
         oceanSilicateFlux, &
         oceanAmmoniumFlux, &
         oceanDMSFlux, &
         oceanDMSPpFlux, &
         oceanDMSPdFlux, &
         oceanHumicsFlux, &
         oceanBlackCarbonFlux, &
         oceanDustIronFlux

    real(kind=RKIND), dimension(:,:), pointer :: &
         oceanAlgaeFlux, &
         oceanDOCFlux, &
         oceanDICFlux, &
         oceanDONFlux, &
         oceanParticulateIronFlux, &
         oceanDissolvedIronFlux

    real(kind=RKIND) :: &
         iceAreaInverse

    integer, pointer :: &
         nCellsSolve, &
         maxAlgaeType, &
         maxDOCType, &
         maxDICType, &
         maxDONType, &
         maxIronType, &
         maxBCType, &
         maxDustType, &
         maxAerosolType

    integer :: &
         iCell, &
         iBioTracers

    call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_zaerosols", config_use_zaerosols)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregate)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolver)
       call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxes)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCoupling)
       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxes)
       call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)
       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

       call MPAS_pool_get_dimension(tracersAggregate, "nCellsSolve", nCellsSolve)

       call MPAS_pool_get_array(tracersAggregate, "iceAreaCell", iceAreaCell)

       call MPAS_pool_get_array(velocitySolver, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocitySolver, "airStressCellV", airStressCellV)

       call MPAS_pool_get_array(atmosFluxes, "sensibleHeatFlux", sensibleHeatFlux)
       call MPAS_pool_get_array(atmosFluxes, "latentHeatFlux", latentHeatFlux)
       call MPAS_pool_get_array(atmosFluxes, "evaporativeWaterFlux", evaporativeWaterFlux)
       call MPAS_pool_get_array(atmosFluxes, "longwaveUp", longwaveUp)

       call MPAS_pool_get_array(shortwave, "absorbedShortwaveFlux", absorbedShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDirectCell", albedoVisibleDirectCell)
       call MPAS_pool_get_array(shortwave, "albedoIRDirectCell", albedoIRDirectCell)
       call MPAS_pool_get_array(shortwave, "albedoVisibleDiffuseCell", albedoVisibleDiffuseCell)
       call MPAS_pool_get_array(shortwave, "albedoIRDiffuseCell", albedoIRDiffuseCell)

       call MPAS_pool_get_array(atmosCoupling, "atmosReferenceHumidity2m", atmosReferenceHumidity2m)
       call MPAS_pool_get_array(atmosCoupling, "atmosReferenceTemperature2m", atmosReferenceTemperature2m)
       call MPAS_pool_get_array(atmosCoupling, "airTemperature", airTemperature)
       call MPAS_pool_get_array(atmosCoupling, "airSpecificHumidity", airSpecificHumidity)

       call MPAS_pool_get_array(oceanCoupling, "seaFreezingTemperature", seaFreezingTemperature)

       call MPAS_pool_get_array(oceanFluxes, "oceanFreshWaterFlux", oceanFreshWaterFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanSaltFlux", oceanSaltFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanHeatFlux", oceanHeatFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanShortwaveFlux", oceanShortwaveFlux)

       call MPAS_pool_get_array(biogeochemistry, "oceanNitrateFlux", oceanNitrateFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanSilicateFlux", oceanSilicateFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanAmmoniumFlux", oceanAmmoniumFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSFlux", oceanDMSFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSPpFlux", oceanDMSPpFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSPdFlux", oceanDMSPdFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanHumicsFlux", oceanHumicsFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDustIronFlux", oceanDustIronFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanBlackCarbonFlux", oceanBlackCarbonFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanAlgaeFlux", oceanAlgaeFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDOCFlux", oceanDOCFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDICFlux", oceanDICFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDONFlux", oceanDONFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanParticulateIronFlux", oceanParticulateIronFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDissolvedIronFlux", oceanDissolvedIronFlux)

       call MPAS_pool_get_dimension(mesh, "maxAlgaeType", maxAlgaeType)
       call MPAS_pool_get_dimension(mesh, "maxDOCType", maxDOCType)
       call MPAS_pool_get_dimension(mesh, "maxDICType", maxDICType)
       call MPAS_pool_get_dimension(mesh, "maxDONType", maxDONType)
       call MPAS_pool_get_dimension(mesh, "maxAerosolType", maxAerosolType)
       call MPAS_pool_get_dimension(mesh, "maxIronType", maxIronType)
       call MPAS_pool_get_dimension(mesh, "maxBCType", maxBCType)
       call MPAS_pool_get_dimension(mesh, "maxDustType", maxDustType)

       do iCell = 1, nCellsSolve

          if (iceAreaCell(iCell) > 0.0_RKIND) then

             iceAreaInverse = 1.0_RKIND / iceAreaCell(iCell)

             airStressCellU(iCell)              = airStressCellU(iCell)              * iceAreaInverse
             airStressCellV(iCell)              = airStressCellV(iCell)              * iceAreaInverse
             sensibleHeatFlux(iCell)            = sensibleHeatFlux(iCell)            * iceAreaInverse
             latentHeatFlux(iCell)              = latentHeatFlux(iCell)              * iceAreaInverse
             absorbedShortwaveFlux(iCell)       = absorbedShortwaveFlux(iCell)       * iceAreaInverse
             longwaveUp(iCell)                  = longwaveUp(iCell)                  * iceAreaInverse
             evaporativeWaterFlux(iCell)        = evaporativeWaterFlux(iCell)        * iceAreaInverse
             atmosReferenceTemperature2m(iCell) = atmosReferenceTemperature2m(iCell) * iceAreaInverse
             atmosReferenceHumidity2m(iCell)    = atmosReferenceHumidity2m(iCell)    * iceAreaInverse
             oceanFreshWaterFlux(iCell)         = oceanFreshWaterFlux(iCell)         * iceAreaInverse
             oceanSaltFlux(iCell)               = oceanSaltFlux(iCell)               * iceAreaInverse
             oceanHeatFlux(iCell)               = oceanHeatFlux(iCell)               * iceAreaInverse
             oceanShortwaveFlux(iCell)          = oceanShortwaveFlux(iCell)          * iceAreaInverse
             albedoVisibleDirectCell(iCell)     = albedoVisibleDirectCell(iCell)     * iceAreaInverse
             albedoIRDirectCell(iCell)          = albedoIRDirectCell(iCell)          * iceAreaInverse
             albedoVisibleDiffuseCell(iCell)    = albedoVisibleDiffuseCell(iCell)    * iceAreaInverse
             albedoIRDiffuseCell(iCell)         = albedoIRDiffuseCell(iCell)         * iceAreaInverse

             if (config_use_zaerosols) then
                oceanDustIronFlux(iCell)       = oceanDustIronFlux(iCell)            * iceAreaInverse
                oceanBlackCarbonFlux(iCell)    = oceanBlackCarbonFlux(iCell)         * iceAreaInverse
             end if

             if (config_use_column_biogeochemistry) then

                oceanNitrateFlux(iCell)        = oceanNitrateFlux(iCell)             * iceAreaInverse
                oceanSilicateFlux(iCell)       = oceanSilicateFlux(iCell)            * iceAreaInverse
                oceanAmmoniumFlux(iCell)       = oceanAmmoniumFlux(iCell)            * iceAreaInverse
                oceanDMSPpFlux(iCell)          = oceanDMSPpFlux(iCell)               * iceAreaInverse
                oceanDMSPdFlux(iCell)          = oceanDMSPdFlux(iCell)               * iceAreaInverse
                oceanDMSFlux(iCell)            = oceanDMSFlux(iCell)                 * iceAreaInverse
                oceanHumicsFlux(iCell)         = oceanHumicsFlux(iCell)              * iceAreaInverse

                do iBioTracers = 1, maxAlgaeType
                   oceanAlgaeFlux(iBioTracers,iCell) = oceanAlgaeFlux(iBioTracers,iCell) * iceAreaInverse
                enddo
                do iBioTracers = 1, maxDOCType
                   oceanDOCFlux(iBioTracers,iCell) = oceanDOCFlux(iBioTracers,iCell) * iceAreaInverse
                enddo
                do iBioTracers = 1, maxDICType
                   oceanDICFlux(iBioTracers,iCell) = oceanDICFlux(iBioTracers,iCell) * iceAreaInverse
                enddo
                do iBioTracers = 1, maxDONType
                   oceanDONFlux(iBioTracers,iCell) = oceanDONFlux(iBioTracers,iCell) * iceAreaInverse
                enddo
                do iBioTracers = 1, maxIronType
                   oceanDissolvedIronFlux(iBioTracers,iCell) = oceanDissolvedIronFlux(iBioTracers,iCell) * iceAreaInverse
                enddo
                do iBioTracers = 1, maxIronType
                   oceanParticulateIronFlux(iBioTracers,iCell) = oceanParticulateIronFlux(iBioTracers,iCell) * iceAreaInverse
                enddo
             endif

          else

             airStressCellU(iCell)              = 0.0_RKIND
             airStressCellV(iCell)              = 0.0_RKIND
             sensibleHeatFlux(iCell)            = 0.0_RKIND
             latentHeatFlux(iCell)              = 0.0_RKIND
             absorbedShortwaveFlux(iCell)       = 0.0_RKIND
             longwaveUp(iCell)                  = &
                  -seaiceStefanBoltzmann * (seaFreezingTemperature(iCell) + seaiceFreshWaterFreezingPoint)**4
             evaporativeWaterFlux(iCell)        = 0.0_RKIND
             atmosReferenceTemperature2m(iCell) = airTemperature(iCell)
             atmosReferenceHumidity2m(iCell)    = airSpecificHumidity(iCell)
             oceanFreshWaterFlux(iCell)         = 0.0_RKIND
             oceanSaltFlux(iCell)               = 0.0_RKIND
             oceanHeatFlux(iCell)               = 0.0_RKIND
             oceanShortwaveFlux(iCell)          = 0.0_RKIND
             albedoVisibleDirectCell(iCell)     = 0.0_RKIND
             albedoIRDirectCell(iCell)          = 0.0_RKIND
             albedoVisibleDiffuseCell(iCell)    = 0.0_RKIND
             albedoIRDiffuseCell(iCell)         = 0.0_RKIND

             if (config_use_zaerosols) then
                oceanDustIronFlux(iCell)        = 0.0_RKIND
                oceanBlackCarbonFlux(iCell)     = 0.0_RKIND
             end if

             if (config_use_column_biogeochemistry) then

                oceanNitrateFlux(iCell)           = 0.0_RKIND
                oceanSilicateFlux(iCell)          = 0.0_RKIND
                oceanAmmoniumFlux(iCell)          = 0.0_RKIND
                oceanDMSPpFlux(iCell)             = 0.0_RKIND
                oceanDMSPdFlux(iCell)             = 0.0_RKIND
                oceanDMSFlux(iCell)               = 0.0_RKIND
                oceanHumicsFlux(iCell)            = 0.0_RKIND
                oceanAlgaeFlux(:,iCell)           = 0.0_RKIND
                oceanDOCFlux(:,iCell)             = 0.0_RKIND
                oceanDICFlux(:,iCell)             = 0.0_RKIND
                oceanDONFlux(:,iCell)             = 0.0_RKIND
                oceanParticulateIronFlux(:,iCell) = 0.0_RKIND
                oceanDissolvedIronFlux(:,iCell)   = 0.0_RKIND

             endif
          endif

       enddo ! iCell

       block => block % next
    enddo

  end subroutine seaice_column_scale_fluxes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack_ocean_mixed_layer
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_icepack_ocean_mixed_layer(domain)

    use icepack_intfc, only: &
         icepack_atm_boundary, &
         icepack_ocn_mixed_layer

    use seaice_constants, only: &
         seaiceOceanAlbedo

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         oceanCoupling, &
         atmosCoupling, &
         atmosForcing, &
         tracersAggregate, &
         drag, &
         oceanFluxes, &
         oceanAtmosphere

    real(kind=RKIND), dimension(:), pointer :: &
         seaSurfaceTemperature, &
         seaFreezingTemperature, &
         oceanMixedLayerDepth, &
         oceanHeatFluxConvergence, &
         airPotentialTemperature, &
         airSpecificHumidity, &
         uAirVelocity, &
         vAirVelocity, &
         windSpeed, &
         airLevelHeight, &
         airDensity, &
         longwaveDown, &
         iceAreaCell, &
         freezingMeltingPotential, &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown, &
         airDragCoefficient, &
         airOceanDragCoefficientRatio, &
         oceanHeatFlux, &
         oceanShortwaveFlux, &
         airStressOceanU, &
         airStressOceanV, &
         atmosReferenceTemperature2mOcean, &
         atmosReferenceHumidity2mOcean, &
         longwaveUpOcean, &
         sensibleHeatFluxOcean, &
         latentHeatFluxOcean, &
         evaporativeWaterFluxOcean, &
         albedoVisibleDirectOcean, &
         albedoIRDirectOcean, &
         albedoVisibleDiffuseOcean, &
         albedoIRDiffuseOcean

    real(kind=RKIND) :: &
         sensibleTransferCoefficient, &
         latentTransferCoefficient, &
         potentialTemperatureDifference, &
         specificHumidityDifference

    real(kind=RKIND), pointer :: &
         config_dt

    integer :: &
         iCell

    integer, pointer :: &
         nCellsSolve

    integer, dimension(:), pointer :: &
           landIceMask

    logical, pointer :: &
         config_use_test_ice_shelf

    call MPAS_pool_get_config(domain % configs, "config_dt", config_dt)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCoupling)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
       call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmosForcing)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregate)
       call MPAS_pool_get_subpool(block % structs, "drag", drag)
       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxes)
       call MPAS_pool_get_subpool(block % structs, "ocean_atmosphere", oceanAtmosphere)

       call MPAS_pool_get_dimension(oceanCoupling, "nCellsSolve", nCellsSolve)

       call MPAS_pool_get_array(oceanCoupling, "seaSurfaceTemperature", seaSurfaceTemperature)
       call MPAS_pool_get_array(oceanCoupling, "seaFreezingTemperature", seaFreezingTemperature)
       call MPAS_pool_get_array(oceanCoupling, "freezingMeltingPotential", freezingMeltingPotential)
       call MPAS_pool_get_array(oceanCoupling, "oceanMixedLayerDepth", oceanMixedLayerDepth)
       call MPAS_pool_get_array(oceanCoupling, "oceanHeatFluxConvergence", oceanHeatFluxConvergence)

       call MPAS_pool_get_array(atmosCoupling, "airPotentialTemperature", airPotentialTemperature)
       call MPAS_pool_get_array(atmosCoupling, "uAirVelocity", uAirVelocity)
       call MPAS_pool_get_array(atmosCoupling, "vAirVelocity", vAirVelocity)
       call MPAS_pool_get_array(atmosCoupling, "airLevelHeight", airLevelHeight)
       call MPAS_pool_get_array(atmosCoupling, "airSpecificHumidity", airSpecificHumidity)
       call MPAS_pool_get_array(atmosCoupling, "airDensity", airDensity)
       call MPAS_pool_get_array(atmosCoupling, "longwaveDown", longwaveDown)
       call MPAS_pool_get_array(atmosCoupling, "shortwaveVisibleDirectDown", shortwaveVisibleDirectDown)
       call MPAS_pool_get_array(atmosCoupling, "shortwaveVisibleDiffuseDown", shortwaveVisibleDiffuseDown)
       call MPAS_pool_get_array(atmosCoupling, "shortwaveIRDirectDown", shortwaveIRDirectDown)
       call MPAS_pool_get_array(atmosCoupling, "shortwaveIRDiffuseDown", shortwaveIRDiffuseDown)

       call MPAS_pool_get_array(atmosForcing, "windSpeed", windSpeed)

       call MPAS_pool_get_array(tracersAggregate, "iceAreaCell", iceAreaCell)

       call MPAS_pool_get_array(drag, "airDragCoefficient", airDragCoefficient)
       call MPAS_pool_get_array(drag, "airOceanDragCoefficientRatio", airOceanDragCoefficientRatio)

       call MPAS_pool_get_array(oceanFluxes, "oceanHeatFlux", oceanHeatFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanShortwaveFlux", oceanShortwaveFlux)

       call MPAS_pool_get_array(oceanAtmosphere, "airStressOceanU", airStressOceanU)
       call MPAS_pool_get_array(oceanAtmosphere, "airStressOceanV", airStressOceanV)
       call MPAS_pool_get_array(oceanAtmosphere, "atmosReferenceTemperature2mOcean", atmosReferenceTemperature2mOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "atmosReferenceHumidity2mOcean", atmosReferenceHumidity2mOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "albedoVisibleDirectOcean", albedoVisibleDirectOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "albedoVisibleDiffuseOcean", albedoVisibleDiffuseOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "albedoIRDirectOcean", albedoIRDirectOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "albedoIRDiffuseOcean", albedoIRDiffuseOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "longwaveUpOcean", longwaveUpOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "sensibleHeatFluxOcean", sensibleHeatFluxOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "latentHeatFluxOcean", latentHeatFluxOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "evaporativeWaterFluxOcean", evaporativeWaterFluxOcean)

       do iCell = 1, nCellsSolve
          call icepack_atm_boundary(&
               sfctype='ocn', &
               Tsf=seaSurfaceTemperature(iCell), &
               potT=airPotentialTemperature(iCell), &
               uatm=uAirVelocity(iCell), &
               vatm=vAirVelocity(iCell), &
               wind=windSpeed(iCell), &
               zlvl=airLevelHeight(iCell), &
               Qa=airSpecificHumidity(iCell), &
               rhoa=airDensity(iCell), &
               strx=airStressOceanU(iCell), &
               stry=airStressOceanV(iCell), &
               Tref=atmosReferenceTemperature2mOcean(iCell), &
               Qref=atmosReferenceHumidity2mOcean(iCell), &
               delt=potentialTemperatureDifference, &
               delq=specificHumidityDifference, &
               lhcoef=latentTransferCoefficient, &
               shcoef=sensibleTransferCoefficient, &
               Cdn_atm=airDragCoefficient(iCell), &
               Cdn_atm_ratio_n=airOceanDragCoefficientRatio(iCell))

          albedoVisibleDirectOcean(iCell)  = seaiceOceanAlbedo
          albedoIRDirectOcean(iCell)       = seaiceOceanAlbedo
          albedoVisibleDiffuseOcean(iCell) = seaiceOceanAlbedo
          albedoIRDiffuseOcean(iCell)      = seaiceOceanAlbedo

          call icepack_ocn_mixed_layer(&
               alvdr_ocn=albedoVisibleDirectOcean(iCell), &
               swvdr=shortwaveVisibleDirectDown(iCell), &
               alidr_ocn=albedoIRDirectOcean(iCell), &
               swidr=shortwaveIRDirectDown(iCell), &
               alvdf_ocn=albedoVisibleDiffuseOcean(iCell), &
               swvdf=shortwaveVisibleDiffuseDown(iCell), &
               alidf_ocn=albedoIRDiffuseOcean(iCell), &
               swidf=shortwaveIRDiffuseDown(iCell), &
               sst=seaSurfaceTemperature(iCell), &
               flwout_ocn=longwaveUpOcean(iCell), &
               fsens_ocn=sensibleHeatFluxOcean(iCell), &
               shcoef=sensibleTransferCoefficient, &
               flat_ocn=latentHeatFluxOcean(iCell), &
               lhcoef=latentTransferCoefficient, &
               evap_ocn=evaporativeWaterFluxOcean(iCell), &
               flw=longwaveDown(iCell), &
               delt=potentialTemperatureDifference, &
               delq=specificHumidityDifference, &
               aice=iceAreaCell(iCell), &
               fhocn=oceanHeatFlux(iCell), &
               fswthru=oceanShortwaveFlux(iCell), &
               hmix=oceanMixedLayerDepth(iCell), &
               Tf=seaFreezingTemperature(iCell), &
               qdp=oceanHeatFluxConvergence(iCell), &
               frzmlt=freezingMeltingPotential(iCell), &
               dt=config_dt)

       enddo ! iCell

       block => block % next
    enddo

    call seaice_icepack_write_warnings(icepack_warnings_aborted())

    ! remove frazil from below ice shelves if we are testing that
    call MPAS_pool_get_config(domain % configs, "config_use_test_ice_shelf", config_use_test_ice_shelf)

    if (config_use_test_ice_shelf) then

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCoupling)

          call MPAS_pool_get_array(oceanCoupling, "freezingMeltingPotential", freezingMeltingPotential)
          call MPAS_pool_get_array(oceanCoupling, "landIceMask", landIceMask)

          call MPAS_pool_get_dimension(oceanCoupling, "nCellsSolve", nCellsSolve)

          do iCell = 1, nCellsSolve

             if (landIceMask(iCell) == 1) then
                freezingMeltingPotential(iCell) = 0.0_RKIND
             endif

          enddo ! iCell

          block => block % next
       enddo

    endif !

  end subroutine seaice_icepack_ocean_mixed_layer

!-----------------------------------------------------------------------
! Other passthrough functions to column package
!-----------------------------------------------------------------------

  subroutine seaice_icepack_init_trcr(&
       airTemperature, &
       seaFreezingTemperature, &
       initialSalinityProfile, &
       initialMeltingTemperatureProfile, &
       surfaceTemperature, &
       iceEnthalpy, &
       snowEnthalpy)

    use icepack_intfc, only: &
         icepack_init_enthalpy

    real(kind=RKIND), intent(in) :: &
         airTemperature, &      ! air temperature (C)
         seaFreezingTemperature ! freezing temperature (C)

    real(kind=RKIND), dimension(:), intent(in) :: &
         initialSalinityProfile, &        ! vertical salinity profile (ppt)
         initialMeltingTemperatureProfile ! vertical temperature profile (C)

    real(kind=RKIND), intent(out) :: &
         surfaceTemperature ! surface temperature (C)

    real(kind=RKIND), dimension(:), intent(out) :: &
         iceEnthalpy, & ! ice enthalpy profile (J/m3)
         snowEnthalpy   ! snow enthalpy profile (J/m3)

    call icepack_init_enthalpy(airTemperature, &
         seaFreezingTemperature, &
         initialSalinityProfile, &
         initialMeltingTemperatureProfile, &
         surfaceTemperature, &
         iceEnthalpy, &
         snowEnthalpy)

    call seaice_icepack_write_warnings(icepack_warnings_aborted())

  end subroutine seaice_icepack_init_trcr

  !-----------------------------------------------------------------------

  subroutine seaice_icepack_init_ocean_conc(&
       oceanAmmoniumConc, &
       oceanDMSPConc, &
       oceanDMSConc, &
       oceanAlgaeConc, &
       oceanDOCConc, &
       oceanDICConc, &
       oceanDONConc, &
       oceanDissolvedIronConc, &
       oceanParticulateIronConc, &
       oceanHumicsConc, &
       oceanNitrateConc, &
       oceanSilicateConc,&
       oceanZAerosolConc, &
       carbonToNitrogenRatioAlgae, &
       carbonToNitrogenRatioDON)

    use icepack_intfc, only: &
         icepack_init_ocean_bio

    real(kind=RKIND), intent(out):: &
         oceanAmmoniumConc, & ! ammonium
         oceanDMSPConc, &     ! DMSPp
         oceanDMSConc, &      ! DMS
         oceanHumicsConc, &   ! humic material
         oceanNitrateConc, &  ! nitrate
         oceanSilicateConc    ! silicate

    real(kind=RKIND), dimension(:), intent(out):: &
         oceanAlgaeConc, &           ! algae
         oceanDOCConc, &             ! DOC
         oceanDICConc, &             ! DIC
         oceanDONConc, &             ! DON
         oceanDissolvedIronConc, &   ! Dissolved Iron
         oceanParticulateIronConc, & ! Particulate Iron
         oceanZAerosolConc           ! BC and dust

    real(kind=RKIND), dimension(:), intent(out), optional :: &
         carbonToNitrogenRatioAlgae, & ! carbon to nitrogen ratio for algae
         carbonToNitrogenRatioDON      ! nitrogen to carbon ratio for proteins

    call icepack_init_ocean_bio(&
         amm=oceanAmmoniumConc, &
         dmsp=oceanDMSPConc, &
         dms=oceanDMSConc, &
         algalN=oceanAlgaeConc(:), &
         doc=oceanDOCConc(:), &
         dic=oceanDICConc(:), &
         don=oceanDONConc(:), &
         fed=oceanDissolvedIronConc(:), &
         fep=oceanParticulateIronConc(:), &
         hum=oceanHumicsConc, &
         nit=oceanNitrateConc, &
         sil=oceanSilicateConc,&
         zaeros=oceanZAerosolConc(:), &
         CToN=carbonToNitrogenRatioAlgae(:), &
         CToN_DON=carbonToNitrogenRatioDON(:))

  end subroutine seaice_icepack_init_ocean_conc

  !-----------------------------------------------------------------------

  function seaice_icepack_enthalpy_ice(iceTemperature, iceSalinity, &
                                       config_thermodynamics_type) result(iceEnthalpy)

    use icepack_intfc, only: icepack_enthalpy_mush

    use seaice_constants, only: &
        seaiceMeltingTemperatureDepression, & ! melting temperature depression factor (C/ppt)
        seaiceDensityIce, &                   ! density of ice (kg/m^3)
        seaiceFreshIceSpecificHeat, &         ! specific heat of fresh ice (J/kg/K)
        seaiceLatentHeatMelting, &            ! latent heat of melting of fresh ice (J/kg)
        seaiceSeaWaterSpecificHeat            ! specific heat of ocn (J/kg/K)

    real(kind=RKIND), intent(in) :: iceTemperature
    real(kind=RKIND), intent(in) :: iceSalinity
    real(kind=RKIND) :: iceEnthalpy

    real(kind=RKIND) :: Tmlt

    character(len=strKIND), pointer :: config_thermodynamics_type

    if (trim(config_thermodynamics_type) == "mushy") then

       iceEnthalpy = icepack_enthalpy_mush(iceTemperature, iceSalinity)

    else

       Tmlt = -iceSalinity*seaiceMeltingTemperatureDepression
       iceEnthalpy = -(seaiceDensityIce * (seaiceFreshIceSpecificHeat*(Tmlt-iceTemperature) &
                     + seaiceLatentHeatMelting*(1.0_RKIND-Tmlt/iceTemperature) &
                     - seaiceSeaWaterSpecificHeat*Tmlt))

    endif

  end function seaice_icepack_enthalpy_ice

!-----------------------------------------------------------------------
! initialize constants
!-----------------------------------------------------------------------

  subroutine seaice_init_icepack_constants()

    use icepack_intfc, only: &
         icepack_configure

    use seaice_constants, only: &
         seaicePuny, &
         iceAreaMinimum, &
         iceThicknessMinimum, &
         snowThicknessMinimum

    call icepack_configure()
    call seaice_icepack_write_warnings(icepack_warnings_aborted())

    iceAreaMinimum       = seaicePuny
    iceThicknessMinimum  = seaicePuny
    snowThicknessMinimum = seaicePuny

  end subroutine seaice_init_icepack_constants

!-----------------------------------------------------------------------
! CICE tracer object
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_tracer_object
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 22nd January 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_tracer_object(domain, tracerObject)

    type(domain_type), intent(in) :: &
         domain

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    integer, pointer :: &
         nCategories, &
         nZBGCTracers

    logical, pointer :: &
         config_use_column_biogeochemistry, &
         config_use_zaerosols

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nCategories", nCategories)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nZBGCTracers", nZBGCTracers)

    call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_zaerosols", config_use_zaerosols)

    ! get the number of CICE tracers in trcrn
    call init_column_tracer_object_tracer_number(domain, tracerObject)

    ! allocate other arrays
    allocate(tracerObject % parentIndex(tracerObject % nTracers))
    allocate(tracerObject % firstAncestorMask(tracerObject % nTracers, tracerObject % nBaseTracers))
    allocate(tracerObject % ancestorIndices(tracerObject % nTracers, tracerObject % nMaxAncestorTracers))
    allocate(tracerObject % ancestorNumber(tracerObject % nTracers))

    ! set the child indices
    call init_column_tracer_object_child_indices(domain, tracerObject)

    ! set the parent indices
    call init_column_tracer_object_parent_indices(domain, tracerObject)

    ! set the first ancestor mask
    call init_column_tracer_object_first_ancestor_mask(domain, tracerObject)

    ! set the ancestor indices
    call init_column_tracer_object_ancestor_indices(domain, tracerObject)

    ! biogeochemistry
    if (config_use_column_biogeochemistry .or. config_use_zaerosols) then

       allocate(tracerObject % index_LayerIndexToDataArray(nZBGCTracers))
       allocate(tracerObject % index_LayerIndexToBioIndex(nZBGCTracers))

       ! set all indices for biogeochemistry including parent, ancestor and ancestor mask
       call init_column_tracer_object_for_biogeochemistry(domain, tracerObject)

    else
       allocate(tracerObject % index_algaeConc(1))
       allocate(tracerObject % index_algalCarbon(1))
       allocate(tracerObject % index_algalChlorophyll(1))
       allocate(tracerObject % index_DOCConc(1))
       allocate(tracerObject % index_DONConc(1))
       allocate(tracerObject % index_DICConc(1))
       allocate(tracerObject % index_dissolvedIronConc(1))
       allocate(tracerObject % index_particulateIronConc(1))
       allocate(tracerObject % index_verticalAerosolsConc(1))

       allocate(tracerObject % index_algaeConcLayer(1))
       allocate(tracerObject % index_algalCarbonLayer(1))
       allocate(tracerObject % index_algalChlorophyllLayer(1))
       allocate(tracerObject % index_DOCConcLayer(1))
       allocate(tracerObject % index_DONConcLayer(1))
       allocate(tracerObject % index_DICConcLayer(1))
       allocate(tracerObject % index_dissolvedIronConcLayer(1))
       allocate(tracerObject % index_particulateIronConcLayer(1))
       allocate(tracerObject % index_verticalAerosolsConcLayer(1))
       allocate(tracerObject % index_verticalAerosolsConcShortwave(1))

       allocate(tracerObject % index_LayerIndexToDataArray(1))
       allocate(tracerObject % index_LayerIndexToBioIndex(1))
    endif

    ! allocate tracer arrays
    !$omp parallel
    allocate(tracerArrayCategory(tracerObject % nTracers, nCategories))
    !$omp end parallel

    allocate(tracerArrayCell(tracerObject % nTracers))

  end subroutine init_column_tracer_object

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_tracer_object_tracer_number
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 22nd January 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_tracer_object_tracer_number(domain, tracerObject)

    type(domain_type), intent(in) :: &
         domain

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols, &
         config_use_brine, &
         config_use_column_biogeochemistry, &
!         config_use_vertical_zsalinity, &      !echmod deprecate
         config_use_vertical_biochemistry, &
         config_use_vertical_tracers, &
         config_use_skeletal_biochemistry, &
         config_use_nitrate, &
         config_use_carbon, &
         config_use_chlorophyll, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_nonreactive, &
         config_use_humics, &
         config_use_DON, &
         config_use_iron, &
         config_use_zaerosols, &
         config_use_effective_snow_density, &
         config_use_snow_grain_radius

    integer, pointer :: &
         nIceLayers, &
         nSnowLayers, &
         nAerosols, &
         nBioLayers, &
         nBioLayersP3, &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON, &
         nParticulateIron, &
         nDissolvedIron, &
         nzAerosols

    integer :: &
         iLayers, &
         iBioTracers, &
         nMobileTracers

    call MPAS_pool_get_config(domain % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(domain % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_aerosols", config_use_aerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_effective_snow_density", config_use_effective_snow_density)
    call MPAS_pool_get_config(domain % configs, "config_use_snow_grain_radius", config_use_snow_grain_radius)

    call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)
!    call MPAS_pool_get_config(domain % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)  !echmod deprecate
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_tracers", config_use_vertical_tracers)
    call MPAS_pool_get_config(domain % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_nitrate", config_use_nitrate)
    call MPAS_pool_get_config(domain % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(domain % configs, "config_use_chlorophyll", config_use_chlorophyll)
    call MPAS_pool_get_config(domain % configs, "config_use_ammonium", config_use_ammonium)
    call MPAS_pool_get_config(domain % configs, "config_use_silicate", config_use_silicate)
    call MPAS_pool_get_config(domain % configs, "config_use_DMS", config_use_DMS)
    call MPAS_pool_get_config(domain % configs, "config_use_nonreactive", config_use_nonreactive)
    call MPAS_pool_get_config(domain % configs, "config_use_humics", config_use_humics)
    call MPAS_pool_get_config(domain % configs, "config_use_DON", config_use_DON)
    call MPAS_pool_get_config(domain % configs, "config_use_iron", config_use_iron)
    call MPAS_pool_get_config(domain % configs, "config_use_zaerosols", config_use_zaerosols)

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nSnowLayers", nSnowLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nAerosols", nAerosols)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nBioLayers", nBioLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nBioLayersP3", nBioLayersP3)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nAlgae", nAlgae)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDOC", nDOC)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDIC", nDIC)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDON", nDON)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nParticulateIron", nParticulateIron)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDissolvedIron", nDissolvedIron)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nzAerosols", nzAerosols)

    !-----------------------------------------------------------------------
    ! physics
    !-----------------------------------------------------------------------

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
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) &
         tracerObject % nTracers = tracerObject % nTracers + 2

    ! level or topo ponds
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) &
         tracerObject % nTracers = tracerObject % nTracers + 1

    ! snow density (density tracer is density from compaction)
    if (config_use_effective_snow_density) then
       tracerObject % nTracers = tracerObject % nTracers + nSnowLayers
    endif

    ! snow grain radius (ice mass, liquid mass, and snow grain radius)
    if (config_use_snow_grain_radius) then
       tracerObject % nTracers = tracerObject % nTracers + nSnowLayers*3
    endif

    ! aerosols
    if (config_use_aerosols) &
         tracerObject % nTracers = tracerObject % nTracers + nAerosols*4

    !-----------------------------------------------------------------------
    ! biogeochemistry
    !-----------------------------------------------------------------------

    if (config_use_column_biogeochemistry .or. config_use_zaerosols) then

       ! save tracer number without bio tracers counted
       tracerObject % nTracersNotBio = tracerObject % nTracers

       ! biogeochemical tracers
       tracerObject % nBioTracersLayer = 0

       ! brine height tracer
       if (config_use_brine) &
            tracerObject % nTracers = tracerObject % nTracers + 1

       ! vertical zSalinity                               !echmod deprecate
!       if (config_use_vertical_zsalinity) then
!          tracerObject % nTracers = tracerObject % nTracers + nBioLayers
!          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + 1
!       endif

       nMobileTracers = 0

       ! Skeletal Biogeochemistry
       if (config_use_skeletal_biochemistry) then
          iLayers = 1
          iBioTracers = 0

          ! Vertical Biogeochemistry
       elseif (config_use_vertical_tracers) then
          iLayers = nBioLayersP3
          iBioTracers = 1

       endif

       ! Algal nitrogen
       if (config_use_vertical_biochemistry .or. config_use_skeletal_biochemistry) then
          tracerObject % nTracers = tracerObject % nTracers + iLayers*nAlgae
          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + nAlgae * iBioTracers
          nMobileTracers = nMobileTracers + nAlgae
       endif

       ! nitrate
       if (config_use_nitrate) then
          tracerObject % nTracers = tracerObject % nTracers + iLayers
          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + 1 * iBioTracers
          nMobileTracers = nMobileTracers + 1
       endif

       ! carbon
       if (config_use_carbon) then
          tracerObject % nTracers = tracerObject % nTracers + iLayers * nDOC &
                                  + iLayers * nDIC
          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + (nDOC + nDIC) * iBioTracers
          nMobileTracers = nMobileTracers + nDIC + nDOC
       endif

       ! Algal chorophyll
       if (config_use_chlorophyll) then
          tracerObject % nTracers = tracerObject % nTracers + iLayers*nAlgae
          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + nAlgae * iBioTracers
          nMobileTracers = nMobileTracers + nAlgae
       endif

       ! ammonium
       if (config_use_ammonium) then
          tracerObject % nTracers = tracerObject % nTracers + iLayers
          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + 1 * iBioTracers
          nMobileTracers = nMobileTracers + 1
       endif

       ! silicate
       if (config_use_silicate) then
          tracerObject % nTracers = tracerObject % nTracers + iLayers
          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + 1 * iBioTracers
          nMobileTracers = nMobileTracers + 1
       endif

       ! DMS
       if (config_use_DMS) then
          tracerObject % nTracers = tracerObject % nTracers + iLayers * 3
          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + 3 * iBioTracers
          nMobileTracers = nMobileTracers + 3
       endif

       ! nonreactive mobile tracer
       if (config_use_nonreactive) then
          tracerObject % nTracers = tracerObject % nTracers + iLayers
          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + 1 * iBioTracers
          nMobileTracers = nMobileTracers + 1
       endif

       ! DON
       if (config_use_DON) then
          tracerObject % nTracers = tracerObject % nTracers + iLayers * nDON
          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + nDON * iBioTracers
          nMobileTracers = nMobileTracers + nDON
       endif

       ! iron
       if (config_use_iron) then
          tracerObject % nTracers = tracerObject % nTracers + iLayers * nParticulateIron &
                                  + iLayers * nDissolvedIron
          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + &
                                            (nParticulateIron + nDissolvedIron) * iBioTracers
          nMobileTracers = nMobileTracers + nParticulateIron + nDissolvedIron
       endif

       ! humic material
       if (config_use_humics) then
          tracerObject % nTracers = tracerObject % nTracers + iLayers
          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + 1 * iBioTracers
          nMobileTracers = nMobileTracers + 1
       endif

       ! zAerosols
       if (config_use_zaerosols) then
          tracerObject % nTracers = tracerObject % nTracers + iLayers * nzAerosols
          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + nzAerosols * iBioTracers
          nMobileTracers = nMobileTracers + nzAerosols
       endif

       ! mobile fraction of vertical tracers
       if (config_use_vertical_tracers) &
            tracerObject % nTracers = tracerObject % nTracers + nMobileTracers

    endif ! config_use_column_biogeochemistry

  end subroutine init_column_tracer_object_tracer_number

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_tracer_object_child_indices
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 22nd January 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_tracer_object_child_indices(domain, tracerObject)

    type(domain_type), intent(in) :: &
         domain

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols, &
         config_use_effective_snow_density, &
         config_use_snow_grain_radius

    integer :: &
         nTracers

    integer, pointer :: &
         nIceLayers, &
         nSnowLayers

    integer, parameter :: indexMissingValue = 0

    call MPAS_pool_get_config(domain % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(domain % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_aerosols", config_use_aerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_effective_snow_density", config_use_effective_snow_density)
    call MPAS_pool_get_config(domain % configs, "config_use_snow_grain_radius", config_use_snow_grain_radius)

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nSnowLayers", nSnowLayers)

    ! ice/snow surface temperature
    tracerObject % index_surfaceTemperature = 1
    nTracers = 1

    ! ice enthalpy
    tracerObject % index_iceEnthalpy = nTracers + 1
    nTracers = nTracers + nIceLayers

    ! snow enthalpy
    tracerObject % index_snowEnthalpy = nTracers + 1
    nTracers = nTracers + nSnowLayers

    ! ice salinity
    tracerObject % index_iceSalinity = nTracers + 1
    nTracers = nTracers + nIceLayers

    ! ice age
    tracerObject % index_iceAge = indexMissingValue
    if (config_use_ice_age) then
       nTracers = nTracers + 1
       tracerObject % index_iceAge = nTracers
    endif

    ! first year ice
    tracerObject % index_firstYearIceArea = indexMissingValue
    if (config_use_first_year_ice) then
       nTracers = nTracers + 1
       tracerObject % index_firstYearIceArea = nTracers
    endif

    ! level ice
    tracerObject % index_levelIceArea   = indexMissingValue
    tracerObject % index_levelIceVolume = indexMissingValue
    if (config_use_level_ice) then
       nTracers = nTracers + 1
       tracerObject % index_levelIceArea = nTracers
       nTracers = nTracers + 1
       tracerObject % index_levelIceVolume = nTracers
    endif

    ! ponds
    tracerObject % index_pondArea         = indexMissingValue
    tracerObject % index_pondDepth        = indexMissingValue
    tracerObject % index_pondLidThickness = indexMissingValue

    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       nTracers = nTracers + 1
       tracerObject % index_pondArea = nTracers
       nTracers = nTracers + 1
       tracerObject % index_pondDepth = nTracers
    endif
    if (config_use_level_meltponds) then
       nTracers = nTracers + 1
       tracerObject % index_pondLidThickness = nTracers
    endif
    if (config_use_topo_meltponds) then
       nTracers = nTracers + 1
       tracerObject % index_pondLidThickness = nTracers
    endif

    ! snow density
    tracerObject % index_snowDensity = indexMissingValue
    if (config_use_effective_snow_density) then
       tracerObject % index_snowDensity = nTracers + 1
       nTracers = nTracers + nSnowLayers
    endif

    ! snow grain radius
    tracerObject % index_snowIceMass = indexMissingValue
    tracerObject % index_snowLiquidMass = indexMissingValue
    tracerObject % index_snowGrainRadius = indexMissingValue
    if (config_use_snow_grain_radius) then
       tracerObject % index_snowIceMass = nTracers + 1
       nTracers = nTracers + nSnowLayers
       tracerObject % index_snowLiquidMass = nTracers + 1
       nTracers = nTracers + nSnowLayers
       tracerObject % index_snowGrainRadius = nTracers + 1
       nTracers = nTracers + nSnowLayers
    endif

    ! aerosols
    tracerObject % index_aerosols = indexMissingValue
    if (config_use_aerosols) then
       tracerObject % index_aerosols = nTracers + 1
    endif

    !-----------------------------------------------------------------------
    ! BGC indices are calculated in the column package
    !-----------------------------------------------------------------------

  end subroutine init_column_tracer_object_child_indices

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_tracer_object_parent_indices
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 22nd January 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_tracer_object_parent_indices(domain, tracerObject)

    type(domain_type), intent(in) :: &
         domain

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols, &
         config_use_effective_snow_density, &
         config_use_snow_grain_radius

    integer :: &
         iIceLayer, &
         iSnowLayer, &
         iAerosol

    integer, pointer :: &
         nIceLayers, &
         nSnowLayers, &
         nAerosols

    call MPAS_pool_get_config(domain % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(domain % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_aerosols", config_use_aerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_effective_snow_density", config_use_effective_snow_density)
    call MPAS_pool_get_config(domain % configs, "config_use_snow_grain_radius", config_use_snow_grain_radius)

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nSnowLayers", nSnowLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nAerosols", nAerosols)

    ! ice/snow surface temperature
    tracerObject % parentIndex(tracerObject % index_surfaceTemperature) = 0

    ! ice enthalpy and salinity
    do iIceLayer = 1, nIceLayers
       tracerObject % parentIndex(tracerObject % index_iceEnthalpy + iIceLayer - 1) = 1
       tracerObject % parentIndex(tracerObject % index_iceSalinity + iIceLayer - 1) = 1
    enddo ! iIceLayer

    ! snow enthalpy
    do iSnowLayer = 1, nSnowLayers
       tracerObject % parentIndex(tracerObject % index_snowEnthalpy + iSnowLayer - 1) = 2
    enddo ! iSnowLayer

    ! ice age
    if (config_use_ice_age) &
         tracerObject % parentIndex(tracerObject % index_iceAge) = 1

    ! first year ice
    if (config_use_first_year_ice) &
         tracerObject % parentIndex(tracerObject % index_firstYearIceArea) = 0

    ! level ice area
    if (config_use_level_ice) then
       tracerObject % parentIndex(tracerObject % index_levelIceArea)   = 0
       tracerObject % parentIndex(tracerObject % index_levelIceVolume) = 1
    endif

    ! level ice ponds
    if (config_use_level_meltponds) then
       tracerObject % parentIndex(tracerObject % index_pondArea)         = 2 + tracerObject % index_levelIceArea
       tracerObject % parentIndex(tracerObject % index_pondDepth)        = 2 + tracerObject % index_pondArea
       tracerObject % parentIndex(tracerObject % index_pondLidThickness) = 2 + tracerObject % index_pondArea
    endif

    ! topo melt ponds
    if (config_use_topo_meltponds) then
       tracerObject % parentIndex(tracerObject % index_pondArea)         = 0
       tracerObject % parentIndex(tracerObject % index_pondDepth)        = 2 + tracerObject % index_pondArea
       tracerObject % parentIndex(tracerObject % index_pondLidThickness) = 2 + tracerObject % index_pondArea
    endif

    ! snow density
    if (config_use_effective_snow_density) then
       do iSnowLayer = 1, nSnowLayers
          tracerObject % parentIndex(tracerObject % index_snowDensity + iSnowLayer - 1) = 2
       enddo ! iSnowLayer
    endif

    ! snow grain radius
    if (config_use_snow_grain_radius) then
       do iSnowLayer = 1, nSnowLayers
          tracerObject % parentIndex(tracerObject % index_snowIceMass + iSnowLayer - 1) = 2
          tracerObject % parentIndex(tracerObject % index_snowLiquidMass + iSnowLayer - 1) = 2
          tracerObject % parentIndex(tracerObject % index_snowGrainRadius + iSnowLayer - 1) = 2
       enddo ! iSnowLayer
    endif

    ! aerosols
    if (config_use_aerosols) then
       do iAerosol = 1, nAerosols
          tracerObject % parentIndex(tracerObject % index_aerosols + (iAerosol-1)*4    ) = 2 ! snow
          tracerObject % parentIndex(tracerObject % index_aerosols + (iAerosol-1)*4 + 1) = 2 ! snow
          tracerObject % parentIndex(tracerObject % index_aerosols + (iAerosol-1)*4 + 2) = 1 ! ice
          tracerObject % parentIndex(tracerObject % index_aerosols + (iAerosol-1)*4 + 3) = 1 ! ice
       enddo ! iAerosol
    endif

    !-----------------------------------------------------------------------
    ! BGC parentIndices are calculated in the column package
    !-----------------------------------------------------------------------

  end subroutine init_column_tracer_object_parent_indices

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_tracer_object_first_ancestor_mask
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 3rd Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_tracer_object_first_ancestor_mask(domain, tracerObject)

    type(domain_type), intent(in) :: &
         domain

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    integer :: &
         iTracer

    ! mask for base quantity on which tracers are carried

    tracerObject % firstAncestorMask = 0.0_RKIND

    do iTracer = 1, tracerObject % nTracers

       if (tracerObject % parentIndex(iTracer) == 0) then

          ! ice area
          tracerObject % firstAncestorMask(iTracer,1) = 1.0_RKIND

       elseif (tracerObject % parentIndex(iTracer) == 1) then  ! ice volume

          ! ice volume
          tracerObject % firstAncestorMask(iTracer,2) = 1.0_RKIND

       elseif (tracerObject % parentIndex(iTracer) == 2) then  ! snow volume

          ! snow volume
          tracerObject % firstAncestorMask(iTracer,3) = 1.0_RKIND

       else

          ! default: ice area
          tracerObject % firstAncestorMask(iTracer,1) = 1.0_RKIND

       endif

    enddo ! iTracer

    !-----------------------------------------------------------------------
    ! BGC firstAncestorMasks are calculated in the column package
    !-----------------------------------------------------------------------

  end subroutine init_column_tracer_object_first_ancestor_mask

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_tracer_object_ancestor_indices
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 3rd Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_tracer_object_ancestor_indices(domain, tracerObject)

    type(domain_type), intent(in) :: &
         domain

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    logical, pointer :: &
         config_use_level_meltponds, &
         config_use_topo_meltponds

    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)

    ! initialize
    tracerObject % ancestorNumber = 0
    tracerObject % ancestorIndices = 0

    ! level melt ponds
    if (config_use_level_meltponds) then

       ! melt pond area
       tracerObject % ancestorNumber (tracerObject % index_pondArea)   = 1
       tracerObject % ancestorIndices(tracerObject % index_pondArea,1) = tracerObject % index_levelIceArea  ! on level ice area

       ! melt pond depth
       tracerObject % ancestorNumber (tracerObject % index_pondDepth)   = 2
       tracerObject % ancestorIndices(tracerObject % index_pondDepth,2) = tracerObject % index_pondArea  ! on melt pond area
       tracerObject % ancestorIndices(tracerObject % index_pondDepth,1) = tracerObject % index_levelIceArea  ! on level ice area

       ! refrozen pond lid
       tracerObject % ancestorNumber (tracerObject % index_pondLidThickness)   = 2
       tracerObject % ancestorIndices(tracerObject % index_pondLidThickness,2) = tracerObject % index_pondArea  ! on melt pond area
       tracerObject % ancestorIndices(tracerObject % index_pondLidThickness,1) = &
            tracerObject % index_levelIceArea  ! on level ice area

    endif

    ! topographic melt ponds
    if (config_use_topo_meltponds) then

       ! melt pond depth
       tracerObject % ancestorNumber (tracerObject % index_pondDepth)   = 1
       tracerObject % ancestorIndices(tracerObject % index_pondDepth,1) = tracerObject % index_pondArea  ! on melt pond area

       ! refrozen pond lid
       tracerObject % ancestorNumber (tracerObject % index_pondLidThickness)   = 1
       tracerObject % ancestorIndices(tracerObject % index_pondLidThickness,1) = tracerObject % index_pondArea  ! on melt pond area

    endif

    !-----------------------------------------------------------------------
    ! BGC ancestorNumbers and ancestorIndices are calculated in the column package
    !-----------------------------------------------------------------------

  end subroutine init_column_tracer_object_ancestor_indices

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  set_cice_tracer_array_category
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 4th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine set_cice_tracer_array_category(block, tracerObject, tracerArrayCategory, iCell, setPhysicsTracers, setBGCTracers)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         tracerArrayCategory

    integer, intent(in) :: &
         iCell

    logical, intent(in) :: &
         setPhysicsTracers, &
         setBGCTracers

    ! get physics tracers
    if (setPhysicsTracers) &
         call set_cice_physics_tracer_array_category(block, tracerArrayCategory, iCell)

    ! get BGC tracers
    if (setBGCTracers) &
         call set_cice_biogeochemistry_tracer_array_category(block, tracerObject, tracerArrayCategory, iCell)

  end subroutine set_cice_tracer_array_category

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_cice_tracer_array_category
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 4th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_cice_tracer_array_category(block, tracerObject, tracerArrayCategory, iCell, getPhysicsTracers, getBGCTracers)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         tracerArrayCategory

    integer, intent(in) :: &
         iCell

    logical, intent(in) :: &
         getPhysicsTracers, &
         getBGCTracers

    ! get physics tracers
    if (getPhysicsTracers) &
         call get_cice_physics_tracer_array_category(block, tracerArrayCategory, iCell)

    ! get BGC tracers
    if (getBGCTracers) &
         call get_cice_biogeochemistry_tracer_array_category(block, tracerObject, tracerArrayCategory, iCell)

  end subroutine get_cice_tracer_array_category

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  set_cice_tracer_array_cell
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 4th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine set_cice_tracer_array_cell(block, tracerObject, tracerArrayCell, iCell, setPhysicsTracers, setBGCTracers)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    real(kind=RKIND), dimension(:), intent(inout) :: &
         tracerArrayCell

    integer, intent(in) :: &
         iCell

    logical, intent(in) :: &
         setPhysicsTracers, &
         setBGCTracers

    ! get physics tracers
    if (setPhysicsTracers) &
         call set_cice_physics_tracer_array_cell(block, tracerArrayCell, iCell)

    ! get BGC tracers
    if (setBGCTracers) &
         call set_cice_biogeochemistry_tracer_array_cell(block, tracerObject, tracerArrayCell, iCell)

  end subroutine set_cice_tracer_array_cell

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_cice_tracer_array_cell
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 4th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_cice_tracer_array_cell(block, tracerObject, tracerArrayCell, iCell, getPhysicsTracers, getBGCTracers)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    real(kind=RKIND), dimension(:), intent(in) :: &
         tracerArrayCell

    integer, intent(in) :: &
         iCell

    logical, intent(in) :: &
         getPhysicsTracers, &
         getBGCTracers

    ! get physics tracers
    if (getPhysicsTracers) &
         call get_cice_physics_tracer_array_cell(block, tracerArrayCell, iCell)

    ! get BGC tracers
    if (getBGCTracers) &
         call get_cice_biogeochemistry_tracer_array_cell(block, tracerObject, tracerArrayCell, iCell)

  end subroutine get_cice_tracer_array_cell

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  set_cice_physics_tracer_array_category
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 4th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine set_cice_physics_tracer_array_category(block, tracerArrayCategory, iCell)

    type(block_type), intent(in) :: &
         block

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         tracerArrayCategory

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols, &
         config_use_effective_snow_density, &
         config_use_snow_grain_radius

    integer, pointer :: &
         nIceLayers, &
         nSnowLayers, &
         nAerosols

    type(MPAS_pool_type), pointer :: &
         tracers

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         surfaceTemperature, &
         iceAge, &
         firstYearIceArea, &
         levelIceArea, &
         levelIceVolume, &
         pondArea, &
         pondDepth, &
         pondLidThickness, &
         iceEnthalpy, &
         snowEnthalpy, &
         iceSalinity, &
         snowScatteringAerosol, &
         snowBodyAerosol, &
         iceScatteringAerosol, &
         iceBodyAerosol, &
         snowIceMass, &
         snowLiquidMass, &
         snowGrainRadius, &
         snowDensity

    integer :: &
         nTracers, &
         iAerosol

    call MPAS_pool_get_config(block % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(block % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(block % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(block % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)
    call MPAS_pool_get_config(block % configs, "config_use_effective_snow_density", config_use_effective_snow_density)
    call MPAS_pool_get_config(block % configs, "config_use_snow_grain_radius", config_use_snow_grain_radius)

    call MPAS_pool_get_dimension(block % dimensions, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(block % dimensions, "nSnowLayers", nSnowLayers)
    call MPAS_pool_get_dimension(block % dimensions, "nAerosols", nAerosols)

    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
    call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity, 1)
    call MPAS_pool_get_array(tracers, "iceAge", iceAge, 1)
    call MPAS_pool_get_array(tracers, "firstYearIceArea", firstYearIceArea, 1)
    call MPAS_pool_get_array(tracers, "levelIceArea", levelIceArea, 1)
    call MPAS_pool_get_array(tracers, "levelIceVolume", levelIceVolume, 1)
    call MPAS_pool_get_array(tracers, "pondArea", pondArea, 1)
    call MPAS_pool_get_array(tracers, "pondDepth", pondDepth, 1)
    call MPAS_pool_get_array(tracers, "pondLidThickness", pondLidThickness, 1)
    call MPAS_pool_get_array(tracers, "snowScatteringAerosol", snowScatteringAerosol, 1)
    call MPAS_pool_get_array(tracers, "snowBodyAerosol", snowBodyAerosol, 1)
    call MPAS_pool_get_array(tracers, "iceScatteringAerosol", iceScatteringAerosol, 1)
    call MPAS_pool_get_array(tracers, "iceBodyAerosol", iceBodyAerosol, 1)
    call MPAS_pool_get_array(tracers, "snowIceMass", snowIceMass, 1)
    call MPAS_pool_get_array(tracers, "snowLiquidMass", snowLiquidMass, 1)
    call MPAS_pool_get_array(tracers, "snowDensity", snowDensity, 1)
    call MPAS_pool_get_array(tracers, "snowGrainRadius", snowGrainRadius, 1)

    nTracers = 1

    ! surfaceTemperature
    tracerArrayCategory(nTracers,:) = surfaceTemperature(1,:,iCell)
    nTracers = nTracers + 1

    ! iceEnthalpy
    tracerArrayCategory(nTracers:nTracers+nIceLayers-1,:) = iceEnthalpy(:,:,iCell)
    nTracers = nTracers + nIceLayers

    ! snowEnthalpy
    tracerArrayCategory(nTracers:nTracers+nSnowLayers-1,:) = snowEnthalpy(:,:,iCell)
    nTracers = nTracers + nSnowLayers

    ! ice Salinity
    tracerArrayCategory(nTracers:nTracers+nIceLayers-1,:) = iceSalinity(:,:,iCell)
    nTracers = nTracers + nIceLayers

    ! iceAge
    if (config_use_ice_age) then
       tracerArrayCategory(nTracers,:) = iceAge(1,:,iCell)
       nTracers = nTracers + 1
    endif

    ! firstYearIceArea
    if (config_use_first_year_ice) then
       tracerArrayCategory(nTracers,:) = firstYearIceArea(1,:,iCell)
       nTracers = nTracers + 1
    endif

    ! level ice tracers
    if (config_use_level_ice) then
       tracerArrayCategory(nTracers,:) = levelIceArea(1,:,iCell)
       nTracers = nTracers + 1
       tracerArrayCategory(nTracers,:) = levelIceVolume(1,:,iCell)
       nTracers = nTracers + 1
    endif

    ! pond tracers
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       tracerArrayCategory(nTracers,:) = pondArea(1,:,iCell)
       nTracers = nTracers + 1
       tracerArrayCategory(nTracers,:) = pondDepth(1,:,iCell)
       nTracers = nTracers + 1
    endif

    ! level or topo ponds
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       tracerArrayCategory(nTracers,:) = pondLidThickness(1,:,iCell)
       nTracers = nTracers + 1
    end if

    ! snow density (density from compaction)
    if (config_use_effective_snow_density) then
       tracerArrayCategory(nTracers:nTracers+nSnowLayers-1,:) = snowDensity(:,:,iCell)
       nTracers = nTracers + nSnowLayers
    endif

    ! snow grain radius
    if (config_use_snow_grain_radius) then
       tracerArrayCategory(nTracers:nTracers+nSnowLayers-1,:) = snowIceMass(:,:,iCell)
       nTracers = nTracers + nSnowLayers
       tracerArrayCategory(nTracers:nTracers+nSnowLayers-1,:) = snowLiquidMass(:,:,iCell)
       nTracers = nTracers + nSnowLayers
       tracerArrayCategory(nTracers:nTracers+nSnowLayers-1,:) = snowGrainRadius(:,:,iCell)
       nTracers = nTracers + nSnowLayers
    endif

    ! aerosols
    if (config_use_aerosols) then
       do iAerosol = 1, nAerosols

          tracerArrayCategory(nTracers+4*(iAerosol-1)  ,:) = snowScatteringAerosol(iAerosol,:,iCell)
          tracerArrayCategory(nTracers+4*(iAerosol-1)+1,:) = snowBodyAerosol(iAerosol,:,iCell)
          tracerArrayCategory(nTracers+4*(iAerosol-1)+2,:) = iceScatteringAerosol(iAerosol,:,iCell)
          tracerArrayCategory(nTracers+4*(iAerosol-1)+3,:) = iceBodyAerosol(iAerosol,:,iCell)

       enddo ! iAerosol
    endif

  end subroutine set_cice_physics_tracer_array_category

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_cice_physics_tracer_array_category
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 4th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_cice_physics_tracer_array_category(block, tracerArrayCategory, iCell)

    type(block_type), intent(inout) :: &
         block

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         tracerArrayCategory

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols, &
         config_use_effective_snow_density, &
         config_use_snow_grain_radius


    integer, pointer :: &
         nIceLayers, &
         nSnowLayers, &
         nAerosols

    type(MPAS_pool_type), pointer :: &
         tracers

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         surfaceTemperature, &
         iceAge, &
         firstYearIceArea, &
         levelIceArea, &
         levelIceVolume, &
         pondArea, &
         pondDepth, &
         pondLidThickness, &
         iceEnthalpy, &
         snowEnthalpy, &
         iceSalinity, &
         snowScatteringAerosol, &
         snowBodyAerosol, &
         iceScatteringAerosol, &
         iceBodyAerosol, &
         snowIceMass, &
         snowLiquidMass, &
         snowGrainRadius, &
         snowDensity

    integer :: &
         nTracers, &
         iAerosol

    call MPAS_pool_get_config(block % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(block % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(block % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(block % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)
    call MPAS_pool_get_config(block % configs, "config_use_effective_snow_density", config_use_effective_snow_density)
    call MPAS_pool_get_config(block % configs, "config_use_snow_grain_radius", config_use_snow_grain_radius)

    call MPAS_pool_get_dimension(block % dimensions, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(block % dimensions, "nSnowLayers", nSnowLayers)
    call MPAS_pool_get_dimension(block % dimensions, "nAerosols", nAerosols)

    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
    call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity, 1)
    call MPAS_pool_get_array(tracers, "iceAge", iceAge, 1)
    call MPAS_pool_get_array(tracers, "firstYearIceArea", firstYearIceArea, 1)
    call MPAS_pool_get_array(tracers, "levelIceArea", levelIceArea, 1)
    call MPAS_pool_get_array(tracers, "levelIceVolume", levelIceVolume, 1)
    call MPAS_pool_get_array(tracers, "pondArea", pondArea, 1)
    call MPAS_pool_get_array(tracers, "pondDepth", pondDepth, 1)
    call MPAS_pool_get_array(tracers, "pondLidThickness", pondLidThickness, 1)
    call MPAS_pool_get_array(tracers, "snowScatteringAerosol", snowScatteringAerosol, 1)
    call MPAS_pool_get_array(tracers, "snowBodyAerosol", snowBodyAerosol, 1)
    call MPAS_pool_get_array(tracers, "iceScatteringAerosol", iceScatteringAerosol, 1)
    call MPAS_pool_get_array(tracers, "iceBodyAerosol", iceBodyAerosol, 1)
    call MPAS_pool_get_array(tracers, "snowIceMass", snowIceMass, 1)
    call MPAS_pool_get_array(tracers, "snowLiquidMass", snowLiquidMass, 1)
    call MPAS_pool_get_array(tracers, "snowDensity", snowDensity, 1)
    call MPAS_pool_get_array(tracers, "snowGrainRadius", snowGrainRadius, 1)

    nTracers = 1

    ! surfaceTemperature
    surfaceTemperature(1,:,iCell) = tracerArrayCategory(nTracers,:)
    nTracers = nTracers + 1

    ! iceEnthalpy
    iceEnthalpy(:,:,iCell) = tracerArrayCategory(nTracers:nTracers+nIceLayers-1,:)
    nTracers = nTracers + nIceLayers

    ! snowEnthalpy
    snowEnthalpy(:,:,iCell) = tracerArrayCategory(nTracers:nTracers+nSnowLayers-1,:)
    nTracers = nTracers + nSnowLayers

    ! ice Salinity
    iceSalinity(:,:,iCell) = tracerArrayCategory(nTracers:nTracers+nIceLayers-1,:)
    nTracers = nTracers + nIceLayers

    ! iceAge
    if (config_use_ice_age) then
       iceAge(1,:,iCell) = tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
    endif

    ! firstYearIceArea
    if (config_use_first_year_ice) then
       firstYearIceArea(1,:,iCell) = tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
    endif

    ! level ice tracers
    if (config_use_level_ice) then
       levelIceArea(1,:,iCell) = tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
       levelIceVolume(1,:,iCell) = tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
    endif

    ! pond tracers
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       pondArea(1,:,iCell) = tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
       pondDepth(1,:,iCell) = tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
    endif

    ! level or topo ponds
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       pondLidThickness(1,:,iCell) = tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
    end if

    ! snow density (density from compaction)
    if (config_use_effective_snow_density) then
       snowDensity(:,:,iCell) = tracerArrayCategory(nTracers:nTracers+nSnowLayers-1,:)
       nTracers = nTracers + nSnowLayers
    endif

    ! snow grain radius
    if (config_use_snow_grain_radius) then
       snowIceMass(:,:,iCell) = tracerArrayCategory(nTracers:nTracers+nSnowLayers-1,:)
       nTracers = nTracers + nSnowLayers
       snowLiquidMass(:,:,iCell) = tracerArrayCategory(nTracers:nTracers+nSnowLayers-1,:)
       nTracers = nTracers + nSnowLayers
       snowGrainRadius(:,:,iCell) = tracerArrayCategory(nTracers:nTracers+nSnowLayers-1,:)
       nTracers = nTracers + nSnowLayers
    endif

    ! aerosols
    if (config_use_aerosols) then
       do iAerosol = 1, nAerosols

          snowScatteringAerosol(iAerosol,:,iCell) = tracerArrayCategory(nTracers+4*(iAerosol-1)  ,:)
          snowBodyAerosol(iAerosol,:,iCell)       = tracerArrayCategory(nTracers+4*(iAerosol-1)+1,:)
          iceScatteringAerosol(iAerosol,:,iCell)  = tracerArrayCategory(nTracers+4*(iAerosol-1)+2,:)
          iceBodyAerosol(iAerosol,:,iCell)        = tracerArrayCategory(nTracers+4*(iAerosol-1)+3,:)

       enddo ! iAerosol
    endif

  end subroutine get_cice_physics_tracer_array_category

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  set_cice_physics_tracer_array_cell
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 4th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine set_cice_physics_tracer_array_cell(block, tracerArrayCell, iCell)

    type(block_type), intent(in) :: &
         block

    real(kind=RKIND), dimension(:), intent(inout) :: &
         tracerArrayCell

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols, &
         config_use_effective_snow_density, &
         config_use_snow_grain_radius

    integer, pointer :: &
         nIceLayers, &
         nSnowLayers, &
         nAerosols

    type(MPAS_pool_type), pointer :: &
         tracers_aggregate

    real(kind=RKIND), dimension(:), pointer :: &
         surfaceTemperatureCell, &
         iceAgeCell, &
         firstYearIceAreaCell, &
         levelIceAreaCell, &
         levelIceVolumeCell, &
         pondAreaCell, &
         pondDepthCell, &
         pondLidThicknessCell

    real(kind=RKIND), dimension(:,:), pointer :: &
         iceEnthalpyCell, &
         snowEnthalpyCell, &
         iceSalinityCell, &
         snowScatteringAerosolCell, &
         snowBodyAerosolCell, &
         iceScatteringAerosolCell, &
         iceBodyAerosolCell, &
         snowIceMassCell, &
         snowLiquidMassCell, &
         snowGrainRadiusCell, &
         snowDensityCell

    integer :: &
         nTracers, &
         iAerosol

    call MPAS_pool_get_config(block % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(block % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(block % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(block % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)
    call MPAS_pool_get_config(block % configs, "config_use_effective_snow_density", config_use_effective_snow_density)
    call MPAS_pool_get_config(block % configs, "config_use_snow_grain_radius", config_use_snow_grain_radius)

    call MPAS_pool_get_dimension(block % dimensions, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(block % dimensions, "nSnowLayers", nSnowLayers)
    call MPAS_pool_get_dimension(block % dimensions, "nAerosols", nAerosols)

    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)

    call MPAS_pool_get_array(tracers_aggregate, "surfaceTemperatureCell", surfaceTemperatureCell)
    call MPAS_pool_get_array(tracers_aggregate, "iceEnthalpyCell", iceEnthalpyCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowEnthalpyCell", snowEnthalpyCell)
    call MPAS_pool_get_array(tracers_aggregate, "iceSalinityCell", iceSalinityCell)
    call MPAS_pool_get_array(tracers_aggregate, "iceAgeCell", iceAgeCell)
    call MPAS_pool_get_array(tracers_aggregate, "firstYearIceAreaCell", firstYearIceAreaCell)
    call MPAS_pool_get_array(tracers_aggregate, "levelIceAreaCell", levelIceAreaCell)
    call MPAS_pool_get_array(tracers_aggregate, "levelIceVolumeCell", levelIceVolumeCell)
    call MPAS_pool_get_array(tracers_aggregate, "pondAreaCell", pondAreaCell)
    call MPAS_pool_get_array(tracers_aggregate, "pondDepthCell", pondDepthCell)
    call MPAS_pool_get_array(tracers_aggregate, "pondLidThicknessCell", pondLidThicknessCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowScatteringAerosolCell", snowScatteringAerosolCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowBodyAerosolCell", snowBodyAerosolCell)
    call MPAS_pool_get_array(tracers_aggregate, "iceScatteringAerosolCell", iceScatteringAerosolCell)
    call MPAS_pool_get_array(tracers_aggregate, "iceBodyAerosolCell", iceBodyAerosolCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowIceMassCell", snowIceMassCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowLiquidMassCell", snowLiquidMassCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowDensityCell", snowDensityCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowGrainRadiusCell", snowGrainRadiusCell)

    nTracers = 1

    ! surfaceTemperature
    tracerArrayCell(nTracers) = surfaceTemperatureCell(iCell)
    nTracers = nTracers + 1

    ! iceEnthalpy
    tracerArrayCell(nTracers:nTracers+nIceLayers-1) = iceEnthalpyCell(:,iCell)
    nTracers = nTracers + nIceLayers

    ! snowEnthalpy
    tracerArrayCell(nTracers:nTracers+nSnowLayers-1) = snowEnthalpyCell(:,iCell)
    nTracers = nTracers + nSnowLayers

    ! ice Salinity
    tracerArrayCell(nTracers:nTracers+nIceLayers-1) = iceSalinityCell(:,iCell)
    nTracers = nTracers + nIceLayers

    ! iceAge
    if (config_use_ice_age) then
       tracerArrayCell(nTracers) = iceAgeCell(iCell)
       nTracers = nTracers + 1
    endif

    ! firstYearIceArea
    if (config_use_first_year_ice) then
       tracerArrayCell(nTracers) = firstYearIceAreaCell(iCell)
       nTracers = nTracers + 1
    endif

    ! level ice tracers
    if (config_use_level_ice) then
       tracerArrayCell(nTracers) = levelIceAreaCell(iCell)
       nTracers = nTracers + 1
       tracerArrayCell(nTracers) = levelIceVolumeCell(iCell)
       nTracers = nTracers + 1
    endif

    ! pond tracers
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       tracerArrayCell(nTracers) = pondAreaCell(iCell)
       nTracers = nTracers + 1
       tracerArrayCell(nTracers) = pondDepthCell(iCell)
       nTracers = nTracers + 1
    endif

    ! level or topo ponds
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       tracerArrayCell(nTracers) = pondLidThicknessCell(iCell)
       nTracers = nTracers + 1
    end if

    ! snow density (density from compaction)
    if (config_use_effective_snow_density) then
       tracerArrayCell(nTracers:nTracers+nSnowLayers-1) = snowDensityCell(:,iCell)
       nTracers = nTracers + nSnowLayers
    endif

    ! snow grain radius
    if (config_use_snow_grain_radius) then
       tracerArrayCell(nTracers:nTracers+nSnowLayers-1) = snowIceMassCell(:,iCell)
       nTracers = nTracers + nSnowLayers
       tracerArrayCell(nTracers:nTracers+nSnowLayers-1) = snowLiquidMassCell(:,iCell)
       nTracers = nTracers + nSnowLayers
       tracerArrayCell(nTracers:nTracers+nSnowLayers-1) = snowGrainRadiusCell(:,iCell)
       nTracers = nTracers + nSnowLayers
    endif

    ! aerosols
    if (config_use_aerosols) then
       do iAerosol = 1, nAerosols

          tracerArrayCell(nTracers+4*(iAerosol-1)  ) = snowScatteringAerosolCell(iAerosol,iCell)
          tracerArrayCell(nTracers+4*(iAerosol-1)+1) = snowBodyAerosolCell(iAerosol,iCell)
          tracerArrayCell(nTracers+4*(iAerosol-1)+2) = iceScatteringAerosolCell(iAerosol,iCell)
          tracerArrayCell(nTracers+4*(iAerosol-1)+3) = iceBodyAerosolCell(iAerosol,iCell)

       enddo ! iAerosol
    endif

  end subroutine set_cice_physics_tracer_array_cell

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_cice_physics_tracer_array_cell
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 4th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_cice_physics_tracer_array_cell(block, tracerArrayCell, iCell)

    type(block_type), intent(inout) :: &
         block

    real(kind=RKIND), dimension(:), intent(in) :: &
         tracerArrayCell

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols, &
         config_use_effective_snow_density, &
         config_use_snow_grain_radius

    integer, pointer :: &
         nIceLayers, &
         nSnowLayers, &
         nAerosols

    type(MPAS_pool_type), pointer :: &
         tracers_aggregate

    real(kind=RKIND), dimension(:), pointer :: &
         surfaceTemperatureCell, &
         iceAgeCell, &
         firstYearIceAreaCell, &
         levelIceAreaCell, &
         levelIceVolumeCell, &
         pondAreaCell, &
         pondDepthCell, &
         pondLidThicknessCell

    real(kind=RKIND), dimension(:,:), pointer :: &
         iceEnthalpyCell, &
         snowEnthalpyCell, &
         iceSalinityCell, &
         snowScatteringAerosolCell, &
         snowBodyAerosolCell, &
         iceScatteringAerosolCell, &
         iceBodyAerosolCell, &
         snowIceMassCell, &
         snowLiquidMassCell, &
         snowGrainRadiusCell, &
         snowDensityCell

    integer :: &
         nTracers, &
         iAerosol

    call MPAS_pool_get_config(block % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(block % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(block % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(block % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)
    call MPAS_pool_get_config(block % configs, "config_use_effective_snow_density", config_use_effective_snow_density)
    call MPAS_pool_get_config(block % configs, "config_use_snow_grain_radius", config_use_snow_grain_radius)

    call MPAS_pool_get_dimension(block % dimensions, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(block % dimensions, "nSnowLayers", nSnowLayers)
    call MPAS_pool_get_dimension(block % dimensions, "nAerosols", nAerosols)

    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)

    call MPAS_pool_get_array(tracers_aggregate, "surfaceTemperatureCell", surfaceTemperatureCell)
    call MPAS_pool_get_array(tracers_aggregate, "iceEnthalpyCell", iceEnthalpyCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowEnthalpyCell", snowEnthalpyCell)
    call MPAS_pool_get_array(tracers_aggregate, "iceSalinityCell", iceSalinityCell)
    call MPAS_pool_get_array(tracers_aggregate, "iceAgeCell", iceAgeCell)
    call MPAS_pool_get_array(tracers_aggregate, "firstYearIceAreaCell", firstYearIceAreaCell)
    call MPAS_pool_get_array(tracers_aggregate, "levelIceAreaCell", levelIceAreaCell)
    call MPAS_pool_get_array(tracers_aggregate, "levelIceVolumeCell", levelIceVolumeCell)
    call MPAS_pool_get_array(tracers_aggregate, "pondAreaCell", pondAreaCell)
    call MPAS_pool_get_array(tracers_aggregate, "pondDepthCell", pondDepthCell)
    call MPAS_pool_get_array(tracers_aggregate, "pondLidThicknessCell", pondLidThicknessCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowScatteringAerosolCell", snowScatteringAerosolCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowBodyAerosolCell", snowBodyAerosolCell)
    call MPAS_pool_get_array(tracers_aggregate, "iceScatteringAerosolCell", iceScatteringAerosolCell)
    call MPAS_pool_get_array(tracers_aggregate, "iceBodyAerosolCell", iceBodyAerosolCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowIceMassCell", snowIceMassCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowLiquidMassCell", snowLiquidMassCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowDensityCell", snowDensityCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowGrainRadiusCell", snowGrainRadiusCell)

    nTracers = 1

    ! surfaceTemperature
    surfaceTemperatureCell(iCell) = tracerArrayCell(nTracers)
    nTracers = nTracers + 1

    ! iceEnthalpy
    iceEnthalpyCell(:,iCell) = tracerArrayCell(nTracers:nTracers+nIceLayers-1)
    nTracers = nTracers + nIceLayers

    ! snowEnthalpy
    snowEnthalpyCell(:,iCell) = tracerArrayCell(nTracers:nTracers+nSnowLayers-1)
    nTracers = nTracers + nSnowLayers

    ! ice Salinity
    iceSalinityCell(:,iCell) = tracerArrayCell(nTracers:nTracers+nIceLayers-1)
    nTracers = nTracers + nIceLayers

    ! iceAge
    if (config_use_ice_age) then
       iceAgeCell(iCell) = tracerArrayCell(nTracers)
       nTracers = nTracers + 1
    endif

    ! firstYearIceArea
    if (config_use_first_year_ice) then
       firstYearIceAreaCell(iCell) = tracerArrayCell(nTracers)
       nTracers = nTracers + 1
    endif

    ! level ice tracers
    if (config_use_level_ice) then
       levelIceAreaCell(iCell) = tracerArrayCell(nTracers)
       nTracers = nTracers + 1
       levelIceVolumeCell(iCell) = tracerArrayCell(nTracers)
       nTracers = nTracers + 1
    endif

    ! pond tracers
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       pondAreaCell(iCell) = tracerArrayCell(nTracers)
       nTracers = nTracers + 1
       pondDepthCell(iCell) = tracerArrayCell(nTracers)
       nTracers = nTracers + 1
    endif

    ! level or topo ponds
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       pondLidThicknessCell(iCell) = tracerArrayCell(nTracers)
       nTracers = nTracers + 1
    end if

    ! snow density (density from compaction)
    if (config_use_effective_snow_density) then
       snowDensityCell(:,iCell) = tracerArrayCell(nTracers:nTracers+nSnowLayers-1)
       nTracers = nTracers + nSnowLayers
    endif

    ! snow grain radius
    if (config_use_snow_grain_radius) then
       snowIceMassCell(:,iCell) = tracerArrayCell(nTracers:nTracers+nSnowLayers-1)
       nTracers = nTracers + nSnowLayers
       snowLiquidMassCell(:,iCell) = tracerArrayCell(nTracers:nTracers+nSnowLayers-1)
       nTracers = nTracers + nSnowLayers
       snowGrainRadiusCell(:,iCell) = tracerArrayCell(nTracers:nTracers+nSnowLayers-1)
       nTracers = nTracers + nSnowLayers
    endif

    ! aerosols
    if (config_use_aerosols) then
       do iAerosol = 1, nAerosols

          snowScatteringAerosolCell(iAerosol,iCell) = tracerArrayCell(nTracers+4*(iAerosol-1)  )
          snowBodyAerosolCell(iAerosol,iCell)       = tracerArrayCell(nTracers+4*(iAerosol-1)+1)
          iceScatteringAerosolCell(iAerosol,iCell)  = tracerArrayCell(nTracers+4*(iAerosol-1)+2)
          iceBodyAerosolCell(iAerosol,iCell)        = tracerArrayCell(nTracers+4*(iAerosol-1)+3)

       enddo ! iAerosol
    endif

  end subroutine get_cice_physics_tracer_array_cell

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  set_cice_biogeochemistry_array_category
!
!> \brief
!> \author Nicole Jeffery
!> \date 12th September 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine set_cice_biogeochemistry_tracer_array_category(block, tracerObject, tracerArrayCategory, iCell)

    type(block_type), intent(in) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         tracerArrayCategory

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_skeletal_biochemistry, &
         config_use_vertical_biochemistry, &
!         config_use_vertical_zsalinity, &      !echmod deprecate
         config_use_vertical_tracers, &
         config_use_brine, &
         config_use_nitrate, &
         config_use_carbon, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_nonreactive, &
         config_use_humics, &
         config_use_DON, &
         config_use_iron, &
         config_use_zaerosols

    integer, pointer :: &
         nBioLayersP3, &
         nBioLayers, &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON, &
         nParticulateIron, &
         nDissolvedIron, &
         nzAerosols

    type(MPAS_pool_type), pointer :: &
         tracers

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         skeletalAlgaeConc, &
         skeletalDOCConc, &
         skeletalDICConc, &
         skeletalDONConc, &
         skeletalDissolvedIronConc, &
         skeletalParticulateIronConc, &
         skeletalNitrateConc, &
         skeletalSilicateConc, &
         skeletalAmmoniumConc, &
         skeletalDMSConc, &
         skeletalDMSPpConc, &
         skeletalDMSPdConc, &
         skeletalNonreactiveConc, &
         skeletalHumicsConc, &
         verticalAlgaeConc, &
         verticalDOCConc, &
         verticalDICConc, &
         verticalDONConc, &
         verticalNitrateConc, &
         verticalSilicateConc, &
         verticalAmmoniumConc, &
         verticalDMSConc, &
         verticalDMSPpConc, &
         verticalDMSPdConc, &
         verticalNonreactiveConc, &
         verticalHumicsConc, &
         verticalParticulateIronConc, &
         verticalDissolvedIronConc, &
         verticalAerosolsConc, &
         verticalSalinity, &
         brineFraction, &
         mobileFraction

    integer :: &
         iBioTracers, &
         iBioCount, &
         iLayers

    call MPAS_pool_get_config(block % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
!    call MPAS_pool_get_config(block % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)      !echmod deprecate
    call MPAS_pool_get_config(block % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_vertical_tracers", config_use_vertical_tracers)
    call MPAS_pool_get_config(block % configs, "config_use_brine", config_use_brine)
    call MPAS_pool_get_config(block % configs, "config_use_nitrate", config_use_nitrate)
    call MPAS_pool_get_config(block % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(block % configs, "config_use_ammonium",config_use_ammonium)
    call MPAS_pool_get_config(block % configs, "config_use_silicate",config_use_silicate)
    call MPAS_pool_get_config(block % configs, "config_use_DMS",config_use_DMS)
    call MPAS_pool_get_config(block % configs, "config_use_nonreactive",config_use_nonreactive)
    call MPAS_pool_get_config(block % configs, "config_use_humics",config_use_humics)
    call MPAS_pool_get_config(block % configs, "config_use_DON",config_use_DON)
    call MPAS_pool_get_config(block % configs, "config_use_iron",config_use_iron)
    call MPAS_pool_get_config(block % configs, "config_use_zaerosols",config_use_zaerosols)

    call MPAS_pool_get_dimension(block % dimensions, "nBioLayers", nBioLayers)
    call MPAS_pool_get_dimension(block % dimensions, "nBioLayersP3", nBioLayersP3)
    call MPAS_pool_get_dimension(block % dimensions, "nzAerosols", nzAerosols)
    call MPAS_pool_get_dimension(block % dimensions, "nAlgae", nAlgae)
    call MPAS_pool_get_dimension(block % dimensions, "nDOC", nDOC)
    call MPAS_pool_get_dimension(block % dimensions, "nDIC", nDIC)
    call MPAS_pool_get_dimension(block % dimensions, "nDON", nDON)
    call MPAS_pool_get_dimension(block % dimensions, "nParticulateIron", nParticulateIron)
    call MPAS_pool_get_dimension(block % dimensions, "nDissolvedIron", nDissolvedIron)

    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_array(tracers, "skeletalAlgaeConc", skeletalAlgaeConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDOCConc", skeletalDOCConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDICConc", skeletalDICConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDONConc", skeletalDONConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalNitrateConc", skeletalNitrateConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalSilicateConc", skeletalSilicateConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalAmmoniumConc", skeletalAmmoniumConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDMSConc", skeletalDMSConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDMSPpConc", skeletalDMSPpConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDMSPdConc", skeletalDMSPdConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalNonreactiveConc", skeletalNonreactiveConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalHumicsConc", skeletalHumicsConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalParticulateIronConc", skeletalParticulateIronConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDissolvedIronConc", skeletalDissolvedIronConc, 1)
    call MPAS_pool_get_array(tracers, "verticalAlgaeConc", verticalAlgaeConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDOCConc", verticalDOCConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDICConc", verticalDICConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDONConc", verticalDONConc, 1)
    call MPAS_pool_get_array(tracers, "verticalNitrateConc", verticalNitrateConc, 1)
    call MPAS_pool_get_array(tracers, "verticalSilicateConc", verticalSilicateConc, 1)
    call MPAS_pool_get_array(tracers, "verticalAmmoniumConc", verticalAmmoniumConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDMSConc", verticalDMSConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDMSPpConc", verticalDMSPpConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDMSPdConc", verticalDMSPdConc, 1)
    call MPAS_pool_get_array(tracers, "verticalNonreactiveConc", verticalNonreactiveConc, 1)
    call MPAS_pool_get_array(tracers, "verticalHumicsConc", verticalHumicsConc, 1)
    call MPAS_pool_get_array(tracers, "verticalParticulateIronConc", verticalParticulateIronConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDissolvedIronConc", verticalDissolvedIronConc, 1)
    call MPAS_pool_get_array(tracers, "verticalAerosolsConc", verticalAerosolsConc, 1)
    call MPAS_pool_get_array(tracers, "verticalSalinity", verticalSalinity, 1)
    call MPAS_pool_get_array(tracers, "brineFraction", brineFraction, 1)
    call MPAS_pool_get_array(tracers, "mobileFraction", mobileFraction, 1)

    ! biogeochemistry

    ! brine height fraction
    if (config_use_brine) &
         tracerArrayCategory(tracerObject % index_brineFraction,:) = brineFraction(1,:,iCell)

    if (config_use_skeletal_biochemistry) then

       ! algal nitrogen
       do iBioTracers = 1, nAlgae
          tracerArrayCategory(tracerObject % index_algaeConc(iBioTracers),:) = &
               skeletalAlgaeConc(iBioTracers,:,iCell)
       enddo

       ! nitrate
       if (config_use_nitrate) &
            tracerArrayCategory(tracerObject % index_nitrateConc,:) = skeletalNitrateConc(1,:,iCell)

       ! DOC
       if (config_use_carbon) then
          do iBioTracers = 1, nDOC
             tracerArrayCategory(tracerObject % index_DOCConc(iBioTracers),:) = skeletalDOCConc(iBioTracers,:,iCell)
          enddo

          ! DIC
          do iBioTracers = 1, nDIC
             tracerArrayCategory(tracerObject % index_DICConc(iBioTracers),:) = skeletalDICConc(iBioTracers,:,iCell)
          enddo
       endif

       ! DON
       if (config_use_DON) then
          do iBioTracers = 1, nDON
             tracerArrayCategory(tracerObject % index_DONConc(iBioTracers),:) = skeletalDONConc(iBioTracers,:,iCell)
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) &
            tracerArrayCategory(tracerObject % index_ammoniumConc,:) = skeletalAmmoniumConc(1,:,iCell)

       ! silicate
       if (config_use_silicate) &
            tracerArrayCategory(tracerObject % index_silicateConc,:) = skeletalSilicateConc(1,:,iCell)

       ! DMS, DMSPp, DMSPd
       if (config_use_DMS) then
          tracerArrayCategory(tracerObject % index_DMSConc,:) = skeletalDMSConc(1,:,iCell)
          tracerArrayCategory(tracerObject % index_DMSPpConc,:) = skeletalDMSPpConc(1,:,iCell)
          tracerArrayCategory(tracerObject % index_DMSPdConc,:) = skeletalDMSPdConc(1,:,iCell)
       endif

       ! nonreactive mobile tracer
       if (config_use_nonreactive) &
            tracerArrayCategory(tracerObject % index_nonreactiveConc,:) = skeletalNonreactiveConc(1,:,iCell)
       ! humic material
       if (config_use_humics) &
            tracerArrayCategory(tracerObject % index_humicsConc,:) = skeletalHumicsConc(1,:,iCell)

       ! Particulate and dissovled Iron
       if (config_use_iron) then
          do iBioTracers = 1, nParticulateIron
             tracerArrayCategory(tracerObject % index_particulateIronConc(iBioTracers),:) = &
                  skeletalParticulateIronConc(iBioTracers,:,iCell)
          enddo
          do iBioTracers = 1, nDissolvedIron
             tracerArrayCategory(tracerObject % index_dissolvedIronConc(iBioTracers),:) = &
                  skeletalDissolvedIronConc(iBioTracers,:,iCell)
          enddo
       endif

    elseif (config_use_vertical_tracers) then

       ! Fraction of biogeochemical tracer in the mobile phase
       do iLayers = 1, tracerObject % nBioTracers
          tracerArrayCategory(tracerObject % index_mobileFraction+iLayers-1,:) = mobileFraction(iLayers,:,iCell)
       enddo

       ! algal nitrogen
       if (config_use_vertical_biochemistry) then
          iBioCount = 0
          do iBioTracers = 1, nAlgae
             do iLayers = 1,nBioLayersP3
                iBiocount = iBiocount + 1
                tracerArrayCategory(tracerObject % index_algaeConc(iBioTracers)+iLayers-1,:) = &
                     verticalAlgaeConc(iBioCount,:,iCell)
             enddo
          enddo
       endif

       ! nitrate
       if (config_use_nitrate) then
          do iLayers = 1, nBioLayersP3
             tracerArrayCategory(tracerObject % index_nitrateConc + iLayers-1,:) = &
                  verticalNitrateConc(iLayers,:,iCell)
          enddo
       endif

       ! DOC
       if (config_use_carbon) then
          iBioCount = 0
          do iBioTracers = 1, nDOC
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerArrayCategory(tracerObject % index_DOCConc(iBioTracers) + iLayers-1,:) = &
                     verticalDOCConc(iBioCount,:,iCell)
             enddo
          enddo
          iBioCount = 0

          ! DIC
          do iBioTracers = 1, nDIC
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerArrayCategory(tracerObject % index_DICConc(iBioTracers) + iLayers-1,:) = &
                     verticalDICConc(iBioCount,:,iCell)
             enddo
          enddo
       endif

       ! DON
       if (config_use_DON) then
          iBioCount = 0
          do iBioTracers = 1, nDON
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerArrayCategory(tracerObject % index_DONConc(iBioTracers) + iLayers-1,:) = &
                     verticalDONConc(iBioCount,:,iCell)
             enddo
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) then
          do iLayers = 1, nBioLayersP3
             tracerArrayCategory(tracerObject % index_ammoniumConc + iLayers-1,:) = &
                  verticalAmmoniumConc(iLayers,:,iCell)
          enddo
       endif

       ! silicate
       if (config_use_silicate) then
          do iLayers = 1, nBioLayersP3
             tracerArrayCategory(tracerObject % index_silicateConc+iLayers-1,:) = &
                  verticalSilicateConc(iLayers,:,iCell)
          enddo
       endif

       ! DMS, DMSPp, DMSPd
       if (config_use_DMS) then
          do iLayers = 1, nBioLayersP3
             tracerArrayCategory(tracerObject % index_DMSConc+iLayers-1,:) = verticalDMSConc(iLayers,:,iCell)
             tracerArrayCategory(tracerObject % index_DMSPpConc+iLayers-1,:) = verticalDMSPpConc(iLayers,:,iCell)
             tracerArrayCategory(tracerObject % index_DMSPdConc+iLayers-1,:) = verticalDMSPdConc(iLayers,:,iCell)
          enddo
       endif

       ! nonreactive purely mobile tracers
       if (config_use_nonreactive) then
          do iLayers = 1, nBioLayersP3
             tracerArrayCategory(tracerObject % index_nonreactiveConc+iLayers-1,:) = &
                  verticalNonreactiveConc(iLayers,:,iCell)
          enddo
       endif

       ! humic material
       if (config_use_humics) then
          do iLayers = 1, nBioLayersP3
             tracerArrayCategory(tracerObject % index_humicsConc+iLayers-1,:) = verticalHumicsConc(iLayers,:,iCell)
          enddo
       endif

       ! particulate and dissolved Iron
       if (config_use_iron) then
          iBioCount = 0
          do iBioTracers = 1, nParticulateIron
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerArrayCategory(tracerObject % index_particulateIronConc(iBioTracers)+iLayers-1,:) = &
                     verticalParticulateIronConc(iBioCount,:,iCell)
             enddo
          enddo
          iBioCount = 0
          do iBioTracers = 1, nDissolvedIron
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerArrayCategory(tracerObject % index_dissolvedIronConc(iBioTracers)+iLayers-1,:) = &
                     verticalDissolvedIronConc(iBioCount,:,iCell)
             enddo
          enddo
       endif

       ! black carbon and dust aerosols
       if (config_use_zaerosols) then
          iBioCount = 0
          do iBioTracers = 1, nzAerosols
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerArrayCategory(tracerObject % index_verticalAerosolsConc(iBioTracers)+iLayers-1,:) = &
                     verticalAerosolsConc(iBioCount,:,iCell)
             enddo
          enddo
       endif

       ! salinity used with BL99 thermodynamics        !echmod deprecate
!       if (config_use_vertical_zsalinity) then
!          do iLayers = 1, nBioLayers
!             tracerArrayCategory(tracerObject % index_verticalSalinity+iLayers-1,:) = &
!                  verticalSalinity(iLayers,:,iCell)
!          enddo
!       endif
    endif

  end subroutine set_cice_biogeochemistry_tracer_array_category

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_cice_biogeochemistry_array_category
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 23rd September 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_cice_biogeochemistry_tracer_array_category(block, tracerObject, tracerArrayCategory, iCell)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         tracerArrayCategory

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_skeletal_biochemistry, &
         config_use_vertical_biochemistry, &
!         config_use_vertical_zsalinity, &      !echmod deprecate
         config_use_vertical_tracers, &
         config_use_brine, &
         config_use_nitrate, &
         config_use_carbon, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_nonreactive, &
         config_use_humics, &
         config_use_DON, &
         config_use_iron, &
         config_use_zaerosols

    integer, pointer :: &
         nBioLayersP3, &
         nBioLayers, &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON, &
         nParticulateIron, &
         nDissolvedIron, &
         nzAerosols

    type(MPAS_pool_type), pointer :: &
         tracers

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         skeletalAlgaeConc, &
         skeletalDOCConc, &
         skeletalDICConc, &
         skeletalDONConc, &
         skeletalDissolvedIronConc, &
         skeletalParticulateIronConc, &
         skeletalNitrateConc, &
         skeletalSilicateConc, &
         skeletalAmmoniumConc, &
         skeletalDMSConc, &
         skeletalDMSPpConc, &
         skeletalDMSPdConc, &
         skeletalNonreactiveConc, &
         skeletalHumicsConc, &
         verticalAlgaeConc, &
         verticalDOCConc, &
         verticalDICConc, &
         verticalDONConc, &
         verticalNitrateConc, &
         verticalSilicateConc, &
         verticalAmmoniumConc, &
         verticalDMSConc, &
         verticalDMSPpConc, &
         verticalDMSPdConc, &
         verticalNonreactiveConc, &
         verticalHumicsConc, &
         verticalParticulateIronConc, &
         verticalDissolvedIronConc, &
         verticalAerosolsConc, &
         verticalSalinity, &
         brineFraction, &
         mobileFraction

    integer :: &
         iBioTracers, &
         iBioCount, &
         iLayers

    call MPAS_pool_get_config(block % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
!    call MPAS_pool_get_config(block % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)       !echmod deprecate
    call MPAS_pool_get_config(block % configs, "config_use_vertical_tracers", config_use_vertical_tracers)
    call MPAS_pool_get_config(block % configs, "config_use_brine", config_use_brine)
    call MPAS_pool_get_config(block % configs, "config_use_nitrate", config_use_nitrate)
    call MPAS_pool_get_config(block % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(block % configs, "config_use_ammonium",config_use_ammonium)
    call MPAS_pool_get_config(block % configs, "config_use_silicate",config_use_silicate)
    call MPAS_pool_get_config(block % configs, "config_use_DMS",config_use_DMS)
    call MPAS_pool_get_config(block % configs, "config_use_nonreactive",config_use_nonreactive)
    call MPAS_pool_get_config(block % configs, "config_use_humics",config_use_humics)
    call MPAS_pool_get_config(block % configs, "config_use_DON",config_use_DON)
    call MPAS_pool_get_config(block % configs, "config_use_iron",config_use_iron)
    call MPAS_pool_get_config(block % configs, "config_use_zaerosols",config_use_zaerosols)

    call MPAS_pool_get_dimension(block % dimensions, "nBioLayers", nBioLayers)
    call MPAS_pool_get_dimension(block % dimensions, "nBioLayersP3", nBioLayersP3)
    call MPAS_pool_get_dimension(block % dimensions, "nzAerosols", nzAerosols)
    call MPAS_pool_get_dimension(block % dimensions, "nAlgae", nAlgae)
    call MPAS_pool_get_dimension(block % dimensions, "nDOC", nDOC)
    call MPAS_pool_get_dimension(block % dimensions, "nDIC", nDIC)
    call MPAS_pool_get_dimension(block % dimensions, "nDON", nDON)
    call MPAS_pool_get_dimension(block % dimensions, "nParticulateIron", nParticulateIron)
    call MPAS_pool_get_dimension(block % dimensions, "nDissolvedIron", nDissolvedIron)

    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_array(tracers, "skeletalAlgaeConc", skeletalAlgaeConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDOCConc", skeletalDOCConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDICConc", skeletalDICConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDONConc", skeletalDONConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalNitrateConc", skeletalNitrateConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalSilicateConc", skeletalSilicateConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalAmmoniumConc", skeletalAmmoniumConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDMSConc", skeletalDMSConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDMSPpConc", skeletalDMSPpConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDMSPdConc", skeletalDMSPdConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalNonreactiveConc", skeletalNonreactiveConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalHumicsConc", skeletalHumicsConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalParticulateIronConc", skeletalParticulateIronConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDissolvedIronConc", skeletalDissolvedIronConc, 1)
    call MPAS_pool_get_array(tracers, "verticalAlgaeConc", verticalAlgaeConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDOCConc", verticalDOCConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDICConc", verticalDICConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDONConc", verticalDONConc, 1)
    call MPAS_pool_get_array(tracers, "verticalNitrateConc", verticalNitrateConc, 1)
    call MPAS_pool_get_array(tracers, "verticalSilicateConc", verticalSilicateConc, 1)
    call MPAS_pool_get_array(tracers, "verticalAmmoniumConc", verticalAmmoniumConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDMSConc", verticalDMSConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDMSPpConc", verticalDMSPpConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDMSPdConc", verticalDMSPdConc, 1)
    call MPAS_pool_get_array(tracers, "verticalNonreactiveConc", verticalNonreactiveConc, 1)
    call MPAS_pool_get_array(tracers, "verticalHumicsConc", verticalHumicsConc, 1)
    call MPAS_pool_get_array(tracers, "verticalParticulateIronConc", verticalParticulateIronConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDissolvedIronConc", verticalDissolvedIronConc, 1)
    call MPAS_pool_get_array(tracers, "verticalAerosolsConc", verticalAerosolsConc, 1)
    call MPAS_pool_get_array(tracers, "verticalSalinity", verticalSalinity, 1)
    call MPAS_pool_get_array(tracers, "brineFraction", brineFraction, 1)
    call MPAS_pool_get_array(tracers, "mobileFraction", mobileFraction, 1)

    ! biogeochemistry
    ! brine height fraction
    if (config_use_brine) &
         brineFraction(1,:,iCell) = tracerArrayCategory(tracerObject % index_brineFraction,:)

    if (config_use_skeletal_biochemistry) then

       ! algal nitrogen
       do iBioTracers = 1, nAlgae
          skeletalAlgaeConc(iBioTracers,:,iCell) = &
               tracerArrayCategory(tracerObject % index_algaeConc(iBioTracers),:)
       enddo

       ! nitrate
       if (config_use_nitrate) &
            skeletalNitrateConc(1,:,iCell) = tracerArrayCategory(tracerObject % index_nitrateConc,:)

       if (config_use_carbon) then

          ! DOC
          do iBioTracers = 1, nDOC
             skeletalDOCConc(iBioTracers,:,iCell) = &
                  tracerArrayCategory(tracerObject % index_DOCConc(iBioTracers),:)
          enddo

          ! DIC
          do iBioTracers = 1, nDIC
             skeletalDICConc(iBioTracers,:,iCell) = &
                  tracerArrayCategory(tracerObject % index_DICConc(iBioTracers),:)
          enddo
       endif

       ! DON
       if (config_use_DON) then
          do iBioTracers = 1, nDON
             skeletalDONConc(iBioTracers,:,iCell) = tracerArrayCategory(tracerObject % index_DONConc(iBioTracers),:)
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) &
            skeletalAmmoniumConc(1,:,iCell) = tracerArrayCategory(tracerObject % index_ammoniumConc,:)
       ! silicate
       if (config_use_silicate) &
            skeletalSilicateConc(1,:,iCell) = tracerArrayCategory(tracerObject % index_silicateConc,:)
       ! DNS, DMSPp, DMSPd
       if (config_use_DMS) then
          skeletalDMSConc(1,:,iCell) = tracerArrayCategory(tracerObject % index_DMSConc,:)
          skeletalDMSPpConc(1,:,iCell) = tracerArrayCategory(tracerObject % index_DMSPpConc,:)
          skeletalDMSPdConc(1,:,iCell) = tracerArrayCategory(tracerObject % index_DMSPdConc,:)
       endif

       ! nonreactive tracer
       if (config_use_nonreactive) &
            skeletalNonreactiveConc(1,:,iCell) = tracerArrayCategory(tracerObject % index_nonreactiveConc,:)
       ! humic material
       if (config_use_humics) &
            skeletalHumicsConc(1,:,iCell) = tracerArrayCategory(tracerObject % index_humicsConc,:)

       if (config_use_iron) then

          ! Particulate Iron
          do iBioTracers = 1, nParticulateIron
             skeletalParticulateIronConc(iBioTracers,:,iCell) = &
                  tracerArrayCategory(tracerObject % index_particulateIronConc(iBioTracers),:)
          enddo

          ! Dissolved Iron
          do iBioTracers = 1, nDissolvedIron
             skeletalDissolvedIronConc(iBioTracers,:,iCell) = &
                  tracerArrayCategory(tracerObject % index_dissolvedIronConc(iBioTracers),:)
          enddo
       endif

    elseif (config_use_vertical_tracers) then

       ! fraction of biogeochemical tracer in the mobile phase
       do iLayers = 1, tracerObject % nBioTracers
          mobileFraction(iLayers,:,iCell) = tracerArrayCategory(tracerObject % index_mobileFraction+iLayers-1,:)
       enddo

       if (config_use_vertical_biochemistry) then
          iBioCount = 0

          ! algal nitrogen
          do iBioTracers = 1, nAlgae
             do iLayers = 1,nBioLayersP3
                iBiocount = iBiocount + 1
                verticalAlgaeConc(iBioCount,:,iCell) = &
                     tracerArrayCategory(tracerObject % index_algaeConc(iBioTracers)+iLayers-1,:)
             enddo
          enddo
       endif

       ! nitrate
       if (config_use_nitrate) then
          do iLayers = 1, nBioLayersP3
             verticalNitrateConc(iLayers,:,iCell) = &
                  tracerArrayCategory(tracerObject % index_nitrateConc + iLayers-1,:)
          enddo
       endif

       if (config_use_carbon) then
          iBioCount = 0

          ! DOC
          do iBioTracers = 1, nDOC
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDOCConc(iBioCount,:,iCell) = &
                     tracerArrayCategory(tracerObject % index_DOCConc(iBioTracers) + iLayers-1,:)
             enddo
          enddo
          iBioCount = 0

          ! DIC
          do iBioTracers = 1, nDIC
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDICConc(iBioCount,:,iCell) = &
                     tracerArrayCategory(tracerObject % index_DICConc(iBioTracers) + iLayers-1,:)
             enddo
          enddo
       endif

       ! DON
       if (config_use_DON) then
          iBioCount = 0
          do iBioTracers = 1, nDON
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDONConc(iBioCount,:,iCell) = &
                     tracerArrayCategory(tracerObject % index_DONConc(iBioTracers) + iLayers-1,:)
             enddo
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) then
          do iLayers = 1, nBioLayersP3
             verticalAmmoniumConc(iLayers,:,iCell) = &
                  tracerArrayCategory(tracerObject % index_ammoniumConc + iLayers-1,:)
          enddo
       endif

       ! silicate
       if (config_use_silicate) then
          do iLayers = 1, nBioLayersP3
             verticalSilicateConc(iLayers,:,iCell) = &
                  tracerArrayCategory(tracerObject % index_silicateConc+iLayers-1,:)
          enddo
       endif

       ! DMS, DMSPp, DMSPd
       if (config_use_DMS) then
          do iLayers = 1, nBioLayersP3
             verticalDMSConc(iLayers,:,iCell) = tracerArrayCategory(tracerObject % index_DMSConc+iLayers-1,:)
             verticalDMSPpConc(iLayers,:,iCell) = tracerArrayCategory(tracerObject % index_DMSPpConc+iLayers-1,:)
             verticalDMSPdConc(iLayers,:,iCell) = tracerArrayCategory(tracerObject % index_DMSPdConc+iLayers-1,:)
          enddo
       endif

       ! nonreactive tracer
       if (config_use_nonreactive) then
          do iLayers = 1, nBioLayersP3
             verticalNonreactiveConc(iLayers,:,iCell) = &
                  tracerArrayCategory(tracerObject % index_nonreactiveConc+iLayers-1,:)
          enddo
       endif

       ! humic material
       if (config_use_humics) then
          do iLayers = 1, nBioLayersP3
             verticalHumicsConc(iLayers,:,iCell) = tracerArrayCategory(tracerObject % index_humicsConc+iLayers-1,:)
          enddo
       endif

       if (config_use_iron) then
          iBioCount = 0

          ! particulate iron
          do iBioTracers = 1, nParticulateIron
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalParticulateIronConc(iBioCount,:,iCell) = &
                     tracerArrayCategory(tracerObject % index_particulateIronConc(iBioTracers)+iLayers-1,:)
             enddo
          enddo
          iBioCount = 0

          ! dissolved iron
          do iBioTracers = 1, nDissolvedIron
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDissolvedIronConc(iBioCount,:,iCell) = &
                     tracerArrayCategory(tracerObject % index_dissolvedIronConc(iBioTracers)+iLayers-1,:)
             enddo
          enddo
       endif

       ! black carbon and dust aerosols
       if (config_use_zaerosols) then
          iBioCount = 0
          do iBioTracers = 1, nzAerosols
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalAerosolsConc(iBioCount,:,iCell) = &
                     tracerArrayCategory(tracerObject % index_verticalAerosolsConc(iBioTracers)+iLayers-1,:)
             enddo
          enddo
       endif

       ! salinity used with BL99 thermodynamics        !echmod deprecate
!       if (config_use_vertical_zsalinity) then
!          do iLayers = 1, nBioLayers
!             verticalSalinity(iLayers,:,iCell) = &
!                  tracerArrayCategory(tracerObject % index_verticalSalinity+iLayers-1,:)
!          enddo
!       endif
    endif

  end subroutine get_cice_biogeochemistry_tracer_array_category

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  set_cice_biogeochemistry_array_cell
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 23rd September 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine set_cice_biogeochemistry_tracer_array_cell(block, tracerObject, tracerArrayCell, iCell)

    type(block_type), intent(in) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    real(kind=RKIND), dimension(:), intent(inout) :: &
         tracerArrayCell

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_skeletal_biochemistry, &
         config_use_vertical_biochemistry, &
!         config_use_vertical_zsalinity, &       !echmod deprecate
         config_use_vertical_tracers, &
         config_use_brine, &
         config_use_nitrate, &
         config_use_carbon, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_nonreactive, &
         config_use_humics, &
         config_use_DON, &
         config_use_iron, &
         config_use_zaerosols

    integer, pointer :: &
         nBioLayersP3, &
         nBioLayersP1, &
         nBioLayers, &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON, &
         nParticulateIron, &
         nDissolvedIron, &
         nzAerosols

    type(MPAS_pool_type), pointer :: &
         tracers_aggregate

    real(kind=RKIND), dimension(:), pointer :: &
         brineFractionCell

    real(kind=RKIND), dimension(:,:), pointer :: &
         skeletalAlgaeConcCell, &
         skeletalDOCConcCell, &
         skeletalDICConcCell, &
         skeletalDONConcCell, &
         skeletalDissolvedIronConcCell, &
         skeletalParticulateIronConcCell, &
         skeletalNitrateConcCell, &
         skeletalSilicateConcCell, &
         skeletalAmmoniumConcCell, &
         skeletalDMSConcCell, &
         skeletalDMSPpConcCell, &
         skeletalDMSPdConcCell, &
         skeletalNonreactiveConcCell, &
         skeletalHumicsConcCell, &
         verticalAlgaeConcCell, &
         verticalDOCConcCell, &
         verticalDICConcCell, &
         verticalDONConcCell, &
         verticalNitrateConcCell, &
         verticalSilicateConcCell, &
         verticalAmmoniumConcCell, &
         verticalDMSConcCell, &
         verticalDMSPpConcCell, &
         verticalDMSPdConcCell, &
         verticalNonreactiveConcCell, &
         verticalHumicsConcCell, &
         verticalParticulateIronConcCell, &
         verticalDissolvedIronConcCell, &
         verticalAerosolsConcCell, &
         verticalSalinityCell, &
         verticalAlgaeIceCell, &
         verticalDOCIceCell, &
         verticalDICIceCell, &
         verticalDONIceCell, &
         verticalNitrateIceCell, &
         verticalSilicateIceCell, &
         verticalAmmoniumIceCell, &
         verticalDMSIceCell, &
         verticalDMSPpIceCell, &
         verticalDMSPdIceCell, &
         verticalNonreactiveIceCell, &
         verticalHumicsIceCell, &
         verticalParticulateIronIceCell, &
         verticalDissolvedIronIceCell, &
         verticalAerosolsIceCell

    integer :: &
         iBioTracers, &
         iBioCount, &
         iLayers, &
         iSnowCount, &
         iIceCount, &
         iBioData

    call MPAS_pool_get_config(block % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
!    call MPAS_pool_get_config(block % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)      !echmod deprecate
    call MPAS_pool_get_config(block % configs, "config_use_vertical_tracers", config_use_vertical_tracers)
    call MPAS_pool_get_config(block % configs, "config_use_brine", config_use_brine)
    call MPAS_pool_get_config(block % configs, "config_use_nitrate", config_use_nitrate)
    call MPAS_pool_get_config(block % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(block % configs, "config_use_ammonium",config_use_ammonium)
    call MPAS_pool_get_config(block % configs, "config_use_silicate",config_use_silicate)
    call MPAS_pool_get_config(block % configs, "config_use_DMS",config_use_DMS)
    call MPAS_pool_get_config(block % configs, "config_use_nonreactive",config_use_nonreactive)
    call MPAS_pool_get_config(block % configs, "config_use_humics",config_use_humics)
    call MPAS_pool_get_config(block % configs, "config_use_DON",config_use_DON)
    call MPAS_pool_get_config(block % configs, "config_use_iron",config_use_iron)
    call MPAS_pool_get_config(block % configs, "config_use_zaerosols",config_use_zaerosols)

    call MPAS_pool_get_dimension(block % dimensions, "nBioLayers", nBioLayers)
    call MPAS_pool_get_dimension(block % dimensions, "nBioLayersP3", nBioLayersP3)
    call MPAS_pool_get_dimension(block % dimensions, "nBioLayersP1", nBioLayersP1)
    call MPAS_pool_get_dimension(block % dimensions, "nzAerosols", nzAerosols)
    call MPAS_pool_get_dimension(block % dimensions, "nAlgae", nAlgae)
    call MPAS_pool_get_dimension(block % dimensions, "nDOC", nDOC)
    call MPAS_pool_get_dimension(block % dimensions, "nDIC", nDIC)
    call MPAS_pool_get_dimension(block % dimensions, "nDON", nDON)
    call MPAS_pool_get_dimension(block % dimensions, "nParticulateIron", nParticulateIron)
    call MPAS_pool_get_dimension(block % dimensions, "nDissolvedIron", nDissolvedIron)

    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)

    call MPAS_pool_get_array(tracers_aggregate, "skeletalAlgaeConcCell", skeletalAlgaeConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDOCConcCell", skeletalDOCConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDICConcCell", skeletalDICConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDONConcCell", skeletalDONConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalNitrateConcCell", skeletalNitrateConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalSilicateConcCell", skeletalSilicateConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalAmmoniumConcCell", skeletalAmmoniumConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDMSConcCell", skeletalDMSConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDMSPpConcCell", skeletalDMSPpConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDMSPdConcCell", skeletalDMSPdConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalNonreactiveConcCell", skeletalNonreactiveConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalHumicsConcCell", skeletalHumicsConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalParticulateIronConcCell", skeletalParticulateIronConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDissolvedIronConcCell", skeletalDissolvedIronConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAlgaeConcCell", verticalAlgaeConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDOCConcCell", verticalDOCConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDICConcCell", verticalDICConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDONConcCell", verticalDONConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalNitrateConcCell", verticalNitrateConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalSilicateConcCell", verticalSilicateConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAmmoniumConcCell", verticalAmmoniumConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDMSConcCell", verticalDMSConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDMSPpConcCell", verticalDMSPpConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDMSPdConcCell", verticalDMSPdConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalNonreactiveConcCell", verticalNonreactiveConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalHumicsConcCell", verticalHumicsConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalParticulateIronConcCell", verticalParticulateIronConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDissolvedIronConcCell", verticalDissolvedIronConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAerosolsConcCell", verticalAerosolsConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAlgaeIceCell", verticalAlgaeIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDOCIceCell", verticalDOCIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDICIceCell", verticalDICIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDONIceCell", verticalDONIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalNitrateIceCell", verticalNitrateIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalSilicateIceCell", verticalSilicateIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAmmoniumIceCell", verticalAmmoniumIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDMSIceCell", verticalDMSIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDMSPpIceCell", verticalDMSPpIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDMSPdIceCell", verticalDMSPdIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalNonreactiveIceCell", verticalNonreactiveIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalHumicsIceCell", verticalHumicsIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalParticulateIronIceCell", verticalParticulateIronIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDissolvedIronIceCell", verticalDissolvedIronIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAerosolsIceCell", verticalAerosolsIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalSalinityCell", verticalSalinityCell)
    call MPAS_pool_get_array(tracers_aggregate, "brineFractionCell", brineFractionCell)

    ! biogeochemistry
    ! brine height fraction
    if (config_use_brine) &
         tracerArrayCell(tracerObject % index_brineFraction)  = brineFractionCell(iCell)

    if (config_use_skeletal_biochemistry) then

       ! algal nitrogen
       do iBioTracers = 1, nAlgae
          tracerArrayCell(tracerObject % index_algaeConc(iBioTracers)) = skeletalAlgaeConcCell(iBioTracers,iCell)
       enddo

       ! nitrate
       if (config_use_nitrate) &
            tracerArrayCell(tracerObject % index_nitrateConc)  = skeletalNitrateConcCell(1,iCell)

       if (config_use_carbon) then

          ! DOC
          do iBioTracers = 1, nDOC
             tracerArrayCell(tracerObject % index_DOCConc(iBioTracers)) = skeletalDOCConcCell(iBioTracers,iCell)
          enddo

          ! DIC
          do iBioTracers = 1, nDIC
             tracerArrayCell(tracerObject % index_DICConc(iBioTracers))  = skeletalDICConcCell(iBioTracers,iCell)
          enddo
       endif

       ! DON
       if (config_use_DON) then
          do iBioTracers = 1, nDON
             tracerArrayCell(tracerObject % index_DONConc(iBioTracers))  = skeletalDONConcCell(iBioTracers,iCell)
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) &
            tracerArrayCell(tracerObject % index_ammoniumConc)  = skeletalAmmoniumConcCell(1,iCell)

       ! silicate
       if (config_use_silicate) &
            tracerArrayCell(tracerObject % index_silicateConc)  = skeletalSilicateConcCell(1,iCell)

       ! DMS, DMSPp, DMSPd
       if (config_use_DMS) then
          tracerArrayCell(tracerObject % index_DMSConc)  = skeletalDMSConcCell(1,iCell)
          tracerArrayCell(tracerObject % index_DMSPpConc) = skeletalDMSPpConcCell(1,iCell)
          tracerArrayCell(tracerObject % index_DMSPdConc) = skeletalDMSPdConcCell(1,iCell)
       endif

       ! nonreactive tracer
       if (config_use_nonreactive) &
            tracerArrayCell(tracerObject % index_nonreactiveConc) = skeletalNonreactiveConcCell(1,iCell)
       ! humic material
       if (config_use_humics) &
            tracerArrayCell(tracerObject % index_humicsConc)  = skeletalHumicsConcCell(1,iCell)

       if (config_use_iron) then

          ! particulate iron
          do iBioTracers = 1, nParticulateIron
             tracerArrayCell(tracerObject % index_particulateIronConc(iBioTracers)) = &
                  skeletalParticulateIronConcCell(iBioTracers,iCell)
          enddo

          ! dissolved iron
          do iBioTracers = 1, nDissolvedIron
             tracerArrayCell(tracerObject % index_dissolvedIronConc(iBioTracers)) = &
                  skeletalDissolvedIronConcCell(iBioTracers,iCell)
          enddo
       endif

    elseif (config_use_vertical_tracers) then

       if (config_use_vertical_biochemistry) then
          iBioCount = 0

          ! algal nitrogen
          do iBioTracers = 1, nAlgae
             iIceCount = (iBioTracers-1)*nBioLayersP1

             do iLayers = 1,nBioLayersP1
                iBiocount = iBiocount + 1
                tracerArrayCell(tracerObject % index_algaeConc(iBioTracers)+iLayers-1) = &
                     verticalAlgaeConcCell(iBioCount,iCell)
                verticalAlgaeIceCell(iLayers+iIceCount,iCell) = verticalAlgaeConcCell(iBioCount,iCell)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBiocount = iBiocount + 1
                tracerArrayCell(tracerObject % index_algaeConc(iBioTracers)+iLayers-1) = &
                     verticalAlgaeConcCell(iBioCount,iCell)
             enddo
          enddo
       endif

       ! nitrate
       if (config_use_nitrate) then
          do iLayers = 1, nBioLayersP1
             tracerArrayCell(tracerObject % index_nitrateConc + iLayers-1) = verticalNitrateConcCell(iLayers,iCell)
             verticalNitrateIceCell(iLayers,iCell) = verticalNitrateConcCell(iLayers,iCell)
          enddo
          do iLayers = nBioLayersP1+1, nBioLayersP3
             tracerArrayCell(tracerObject % index_nitrateConc + iLayers-1) = verticalNitrateConcCell(iLayers,iCell)
          enddo
       endif

       if (config_use_carbon) then
          iBioCount = 0

          ! DOC
          do iBioTracers = 1, nDOC
             iIceCount = (iBioTracers-1)*nBioLayersP1

             do iLayers = 1,nBioLayersP1
                iBioCount = iBioCount + 1
                tracerArrayCell(tracerObject % index_DOCConc(iBioTracers) + iLayers-1) = &
                     verticalDOCConcCell(iBioCount,iCell)
                verticalDOCIceCell(iLayers+iIceCount,iCell) = verticalDOCConcCell(iBioCount,iCell)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerArrayCell(tracerObject % index_DOCConc(iBioTracers) + iLayers-1) = &
                     verticalDOCConcCell(iBioCount,iCell)
             enddo
          enddo
          iBioCount = 0

          ! DIC
          do iBioTracers = 1, nDIC
             iIceCount = (iBioTracers-1)*nBioLayersP1

             do iLayers = 1,nBioLayersP1
                iBioCount = iBioCount + 1
                tracerArrayCell(tracerObject % index_DICConc(iBioTracers) + iLayers-1) = &
                     verticalDICConcCell(iBioCount,iCell)
                verticalDICIceCell(iLayers+iIceCount,iCell) = verticalDICConcCell(iBioCount,iCell)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerArrayCell(tracerObject % index_DICConc(iBioTracers) + iLayers-1) = &
                     verticalDICConcCell(iBioCount,iCell)
             enddo
          enddo
       endif

       ! DON
       if (config_use_DON) then
          iBioCount = 0
          do iBioTracers = 1, nDON
             iIceCount = (iBioTracers-1)*nBioLayersP1

             do iLayers = 1,nBioLayersP1
                iBioCount = iBioCount + 1
                tracerArrayCell(tracerObject % index_DONConc(iBioTracers) + iLayers-1) = &
                     verticalDONConcCell(iBioCount,iCell)
                verticalDONIceCell(iLayers+iIceCount,iCell) = verticalDONConcCell(iBioCount,iCell)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerArrayCell(tracerObject % index_DONConc(iBioTracers) + iLayers-1) = &
                     verticalDONConcCell(iBioCount,iCell)
             enddo
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) then
          do iLayers = 1, nBioLayersP1
             tracerArrayCell(tracerObject % index_ammoniumConc + iLayers-1) = &
                  verticalAmmoniumConcCell(iLayers,iCell)
             verticalAmmoniumIceCell(iLayers,iCell) = verticalAmmoniumConcCell(iLayers,iCell)
          enddo
          do iLayers = nBioLayersP1+1, nBioLayersP3
             tracerArrayCell(tracerObject % index_ammoniumConc + iLayers-1) = &
                  verticalAmmoniumConcCell(iLayers,iCell)
          enddo
       endif

       ! silicate
       if (config_use_silicate) then
          do iLayers = 1, nBioLayersP1
             tracerArrayCell(tracerObject % index_silicateConc+iLayers-1) = verticalSilicateConcCell(iLayers,iCell)
             verticalSilicateIceCell(iLayers,iCell) = verticalSilicateConcCell(iLayers,iCell)
          enddo
          do iLayers = nBioLayersP1+1, nBioLayersP3
             tracerArrayCell(tracerObject % index_silicateConc+iLayers-1) = verticalSilicateConcCell(iLayers,iCell)
          enddo
       endif

       ! DMS, DMSPp, DMSPd
       if (config_use_DMS) then
          do iLayers = 1, nBioLayersP1
             tracerArrayCell(tracerObject % index_DMSConc+iLayers-1)  = verticalDMSConcCell(iLayers,iCell)
             tracerArrayCell(tracerObject % index_DMSPpConc+iLayers-1)  = verticalDMSPpConcCell(iLayers,iCell)
             tracerArrayCell(tracerObject % index_DMSPdConc+iLayers-1)  = verticalDMSPdConcCell(iLayers,iCell)
             verticalDMSIceCell(iLayers,iCell)  = verticalDMSConcCell(iLayers,iCell)
             verticalDMSPpIceCell(iLayers,iCell)  = verticalDMSPpConcCell(iLayers,iCell)
             verticalDMSPdIceCell(iLayers,iCell)  = verticalDMSPdConcCell(iLayers,iCell)
          enddo
          do iLayers = nBioLayersP1+1, nBioLayersP3
             tracerArrayCell(tracerObject % index_DMSConc+iLayers-1)  = verticalDMSConcCell(iLayers,iCell)
             tracerArrayCell(tracerObject % index_DMSPpConc+iLayers-1)  = verticalDMSPpConcCell(iLayers,iCell)
             tracerArrayCell(tracerObject % index_DMSPdConc+iLayers-1)  = verticalDMSPdConcCell(iLayers,iCell)
          enddo
       endif

       ! nonreactive
       if (config_use_nonreactive) then
          do iLayers = 1, nBioLayersP1
             tracerArrayCell(tracerObject % index_nonreactiveConc+iLayers-1) = &
                  verticalNonreactiveConcCell(iLayers,iCell)
             verticalNonreactiveIceCell(iLayers,iCell) = verticalNonreactiveConcCell(iLayers,iCell)
          enddo
          do iLayers = nBioLayersP1+1, nBioLayersP3
             tracerArrayCell(tracerObject % index_nonreactiveConc+iLayers-1) = &
                  verticalNonreactiveConcCell(iLayers,iCell)
          enddo
       endif

       ! humic material
       if (config_use_humics) then
          do iLayers = 1, nBioLayersP1
             tracerArrayCell(tracerObject % index_humicsConc+iLayers-1)  = verticalHumicsConcCell(iLayers,iCell)
             verticalHumicsIceCell(iLayers,iCell) = verticalHumicsConcCell(iLayers,iCell)
          enddo
          do iLayers = nBioLayersP1+1, nBioLayersP3
             tracerArrayCell(tracerObject % index_humicsConc+iLayers-1)  = verticalHumicsConcCell(iLayers,iCell)
          enddo
       endif

       if (config_use_iron) then
          iBioCount = 0

          ! particulate iron
          do iBioTracers = 1, nParticulateIron
             iIceCount = (iBioTracers-1)*nBioLayersP1

             do iLayers = 1,nBioLayersP1
                iBioCount = iBioCount + 1
                tracerArrayCell(tracerObject % index_particulateIronConc(iBioTracers)+iLayers-1) = &
                     verticalParticulateIronConcCell(iBioCount,iCell)
                verticalParticulateIronIceCell(iLayers+iIceCount,iCell) = verticalParticulateIronConcCell(iBioCount,iCell)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerArrayCell(tracerObject % index_particulateIronConc(iBioTracers)+iLayers-1) = &
                     verticalParticulateIronConcCell(iBioCount,iCell)
             enddo
          enddo
          iBioCount = 0

          ! dissolved iron
          do iBioTracers = 1, nDissolvedIron
             iIceCount = (iBioTracers-1)*nBioLayersP1

             do iLayers = 1,nBioLayersP1
                iBioCount = iBioCount + 1
                tracerArrayCell(tracerObject % index_dissolvedIronConc(iBioTracers)+iLayers-1) = &
                     verticalDissolvedIronConcCell(iBioCount,iCell)
                verticalDissolvedIronIceCell(iLayers+iIceCount,iCell) = verticalDissolvedIronConcCell(iBioCount,iCell)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerArrayCell(tracerObject % index_dissolvedIronConc(iBioTracers)+iLayers-1) = &
                     verticalDissolvedIronConcCell(iBioCount,iCell)
             enddo
          enddo
       endif

       ! black carbon and dust aerosols
       if (config_use_zaerosols) then
          iBioCount = 0
          do iBioTracers = 1, nzAerosols
             iIceCount = (iBioTracers-1)*nBioLayersP1

             do iLayers = 1,nBioLayersP1
                iBioCount = iBioCount + 1
                tracerArrayCell(tracerObject % index_verticalAerosolsConc(iBioTracers)+iLayers-1) = &
                     verticalAerosolsConcCell(iBioCount,iCell)
                verticalAerosolsIceCell(iLayers+iIceCount,iCell) = verticalAerosolsConcCell(iBioCount,iCell)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerArrayCell(tracerObject % index_verticalAerosolsConc(iBioTracers)+iLayers-1) = &
                     verticalAerosolsConcCell(iBioCount,iCell)
             enddo
          enddo
       endif

       ! salinity used with BL99 thermodynamics        !echmod deprecate
!       if (config_use_vertical_zsalinity) then
!          do iLayers = 1, nBioLayers
!             tracerArrayCell(tracerObject % index_verticalSalinity+iLayers-1)  = verticalSalinityCell(iLayers,iCell)
!          enddo
!       endif
    endif

  end subroutine set_cice_biogeochemistry_tracer_array_cell

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_cice_biogeochemistry_tracer_array_cell
!
!> \brief
!> \author Nicole Jeffery
!> \date 23rd September 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_cice_biogeochemistry_tracer_array_cell(block, tracerObject, tracerArrayCell, iCell)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    real(kind=RKIND), dimension(:), intent(in) :: &
         tracerArrayCell

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_skeletal_biochemistry, &
         config_use_vertical_biochemistry, &
!         config_use_vertical_zsalinity, &       !echmod deprecate
         config_use_vertical_tracers, &
         config_use_brine, &
         config_use_nitrate, &
         config_use_carbon, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_nonreactive, &
         config_use_humics, &
         config_use_DON, &
         config_use_iron, &
         config_use_zaerosols

    integer, pointer :: &
         nBioLayersP3, &
         nBioLayersP1, &
         nBioLayers, &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON, &
         nParticulateIron, &
         nDissolvedIron, &
         nzAerosols

    type(MPAS_pool_type), pointer :: &
         tracers_aggregate, &
         biogeochemistry

    real(kind=RKIND), dimension(:), pointer :: &
         brineFractionCell, &
         carbonToNitrogenRatioAlgae, &
         carbonToNitrogenRatioDON

    real(kind=RKIND), dimension(:,:), pointer :: &
         skeletalAlgaeConcCell, &
         skeletalDOCConcCell, &
         skeletalDICConcCell, &
         skeletalDONConcCell, &
         skeletalDissolvedIronConcCell, &
         skeletalParticulateIronConcCell, &
         skeletalNitrateConcCell, &
         skeletalSilicateConcCell, &
         skeletalAmmoniumConcCell, &
         skeletalDMSConcCell, &
         skeletalDMSPpConcCell, &
         skeletalDMSPdConcCell, &
         skeletalNonreactiveConcCell, &
         skeletalHumicsConcCell, &
         verticalAlgaeConcCell, &
         verticalDOCConcCell, &
         verticalDICConcCell, &
         verticalDONConcCell, &
         verticalNitrateConcCell, &
         verticalSilicateConcCell, &
         verticalAmmoniumConcCell, &
         verticalDMSConcCell, &
         verticalDMSPpConcCell, &
         verticalDMSPdConcCell, &
         verticalNonreactiveConcCell, &
         verticalHumicsConcCell, &
         verticalParticulateIronConcCell, &
         verticalDissolvedIronConcCell, &
         verticalAerosolsConcCell, &
         verticalSalinityCell, &
         verticalAlgaeIceCell, &
         verticalDOCIceCell, &
         verticalDICIceCell, &
         verticalDONIceCell, &
         verticalNitrateIceCell, &
         verticalSilicateIceCell, &
         verticalAmmoniumIceCell, &
         verticalDMSIceCell, &
         verticalDMSPpIceCell, &
         verticalDMSPdIceCell, &
         verticalNonreactiveIceCell, &
         verticalHumicsIceCell, &
         verticalParticulateIronIceCell, &
         verticalDissolvedIronIceCell, &
         verticalAerosolsIceCell, &
         verticalAerosolsSnowCell, &
         verticalDOCLabileIceCell, &
         verticalAlgaeTotalCarbonIceCell, &
         verticalBCTotalIceCell, &
         verticalDustTotalIceCell, &
         verticalBCTotalSnowCell, &
         verticalDustTotalSnowCell

    integer :: &
         iBioTracers, &
         iBioCount, &
         iLayers, &
         iIceCount, &
         iSnowCount, &
         nAeroType

    call MPAS_pool_get_config(block % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
!    call MPAS_pool_get_config(block % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)       !echmod deprecate
    call MPAS_pool_get_config(block % configs, "config_use_vertical_tracers", config_use_vertical_tracers)
    call MPAS_pool_get_config(block % configs, "config_use_brine", config_use_brine)
    call MPAS_pool_get_config(block % configs, "config_use_nitrate", config_use_nitrate)
    call MPAS_pool_get_config(block % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(block % configs, "config_use_ammonium",config_use_ammonium)
    call MPAS_pool_get_config(block % configs, "config_use_silicate",config_use_silicate)
    call MPAS_pool_get_config(block % configs, "config_use_DMS",config_use_DMS)
    call MPAS_pool_get_config(block % configs, "config_use_nonreactive",config_use_nonreactive)
    call MPAS_pool_get_config(block % configs, "config_use_humics",config_use_humics)
    call MPAS_pool_get_config(block % configs, "config_use_DON",config_use_DON)
    call MPAS_pool_get_config(block % configs, "config_use_iron",config_use_iron)
    call MPAS_pool_get_config(block % configs, "config_use_zaerosols",config_use_zaerosols)

    call MPAS_pool_get_dimension(block % dimensions, "nBioLayers", nBioLayers)
    call MPAS_pool_get_dimension(block % dimensions, "nBioLayersP3", nBioLayersP3)
    call MPAS_pool_get_dimension(block % dimensions, "nBioLayersP1", nBioLayersP1)
    call MPAS_pool_get_dimension(block % dimensions, "nzAerosols", nzAerosols)
    call MPAS_pool_get_dimension(block % dimensions, "nAlgae", nAlgae)
    call MPAS_pool_get_dimension(block % dimensions, "nDOC", nDOC)
    call MPAS_pool_get_dimension(block % dimensions, "nDIC", nDIC)
    call MPAS_pool_get_dimension(block % dimensions, "nDON", nDON)
    call MPAS_pool_get_dimension(block % dimensions, "nParticulateIron", nParticulateIron)
    call MPAS_pool_get_dimension(block % dimensions, "nDissolvedIron", nDissolvedIron)

    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
    call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)

    call MPAS_pool_get_array(tracers_aggregate, "skeletalAlgaeConcCell", skeletalAlgaeConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDOCConcCell", skeletalDOCConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDICConcCell", skeletalDICConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDONConcCell", skeletalDONConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalNitrateConcCell", skeletalNitrateConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalSilicateConcCell", skeletalSilicateConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalAmmoniumConcCell", skeletalAmmoniumConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDMSConcCell", skeletalDMSConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDMSPpConcCell", skeletalDMSPpConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDMSPdConcCell", skeletalDMSPdConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalNonreactiveConcCell", skeletalNonreactiveConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalHumicsConcCell", skeletalHumicsConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalParticulateIronConcCell", skeletalParticulateIronConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDissolvedIronConcCell", skeletalDissolvedIronConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAlgaeConcCell", verticalAlgaeConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDOCConcCell", verticalDOCConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDICConcCell", verticalDICConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDONConcCell", verticalDONConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalNitrateConcCell", verticalNitrateConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalSilicateConcCell", verticalSilicateConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAmmoniumConcCell", verticalAmmoniumConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDMSConcCell", verticalDMSConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDMSPpConcCell", verticalDMSPpConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDMSPdConcCell", verticalDMSPdConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalNonreactiveConcCell", verticalNonreactiveConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalHumicsConcCell", verticalHumicsConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalParticulateIronConcCell", verticalParticulateIronConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDissolvedIronConcCell", verticalDissolvedIronConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAerosolsConcCell", verticalAerosolsConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAlgaeIceCell", verticalAlgaeIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDOCIceCell", verticalDOCIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDICIceCell", verticalDICIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDONIceCell", verticalDONIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalNitrateIceCell", verticalNitrateIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalSilicateIceCell", verticalSilicateIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAmmoniumIceCell", verticalAmmoniumIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDMSIceCell", verticalDMSIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDMSPpIceCell", verticalDMSPpIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDMSPdIceCell", verticalDMSPdIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalNonreactiveIceCell", verticalNonreactiveIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalHumicsIceCell", verticalHumicsIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalParticulateIronIceCell", verticalParticulateIronIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDissolvedIronIceCell", verticalDissolvedIronIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAerosolsIceCell", verticalAerosolsIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAerosolsSnowCell", verticalAerosolsSnowCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalSalinityCell", verticalSalinityCell)
    call MPAS_pool_get_array(tracers_aggregate, "brineFractionCell", brineFractionCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAlgaeTotalCarbonIceCell", verticalAlgaeTotalCarbonIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDOCLabileIceCell", verticalDOCLabileIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalBCTotalIceCell", verticalBCTotalIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDustTotalIceCell", verticalDustTotalIceCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDustTotalSnowCell", verticalDustTotalSnowCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalBCTotalSnowCell", verticalBCTotalSnowCell)
    call MPAS_pool_get_array(biogeochemistry, "carbonToNitrogenRatioAlgae", carbonToNitrogenRatioAlgae)
    call MPAS_pool_get_array(biogeochemistry, "carbonToNitrogenRatioDON", carbonToNitrogenRatioDON)

    ! biogeochemistry
    ! brine height fraction
    if (config_use_brine) &
         brineFractionCell(iCell) = tracerArrayCell(tracerObject % index_brineFraction)

    if (config_use_skeletal_biochemistry) then

       ! algal nitrogen
       do iBioTracers = 1, nAlgae
          skeletalAlgaeConcCell(iBioTracers,iCell) = tracerArrayCell(tracerObject % index_algaeConc(iBioTracers))
       enddo

       ! nitrate
       if (config_use_nitrate) &
            skeletalNitrateConcCell(1,iCell) = tracerArrayCell(tracerObject % index_nitrateConc)

       if (config_use_carbon) then

          ! DOC
          do iBioTracers = 1, nDOC
             skeletalDOCConcCell(iBioTracers,iCell) =  tracerArrayCell(tracerObject % index_DOCConc(iBioTracers))
          enddo

          ! DIC
          do iBioTracers = 1, nDIC
             skeletalDICConcCell(iBioTracers,iCell) = tracerArrayCell(tracerObject % index_DICConc(iBioTracers))
          enddo
       endif

       ! DON
       if (config_use_DON) then
          do iBioTracers = 1, nDON
             skeletalDONConcCell(iBioTracers,iCell) = tracerArrayCell(tracerObject % index_DONConc(iBioTracers))
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) &
            skeletalAmmoniumConcCell(1,iCell) = tracerArrayCell(tracerObject % index_ammoniumConc)

       ! silicate
       if (config_use_silicate) &
            skeletalSilicateConcCell(1,iCell) = tracerArrayCell(tracerObject % index_silicateConc)

       ! DMS
       if (config_use_DMS) then
          skeletalDMSConcCell(1,iCell) = tracerArrayCell(tracerObject % index_DMSConc)
          skeletalDMSPpConcCell(1,iCell) = tracerArrayCell(tracerObject % index_DMSPpConc)
          skeletalDMSPdConcCell(1,iCell) = tracerArrayCell(tracerObject % index_DMSPdConc)
       endif

       ! nonreactive tracer
       if (config_use_nonreactive) &
            skeletalNonreactiveConcCell(1,iCell) = tracerArrayCell(tracerObject % index_nonreactiveConc)
       ! humic material
       if (config_use_humics) &
            skeletalHumicsConcCell(1,iCell) = tracerArrayCell(tracerObject % index_humicsConc)

       if (config_use_iron) then

          ! particulate iron
          do iBioTracers = 1, nParticulateIron
             skeletalParticulateIronConcCell(iBioTracers,iCell) = &
                  tracerArrayCell(tracerObject % index_particulateIronConc(iBioTracers))
          enddo

          ! dissolved iron
          do iBioTracers = 1, nDissolvedIron
             skeletalDissolvedIronConcCell(iBioTracers,iCell) = &
                  tracerArrayCell(tracerObject % index_dissolvedIronConc(iBioTracers))
          enddo
       endif

    elseif (config_use_vertical_tracers) then

       if (config_use_vertical_biochemistry) then
          iBioCount = 0

          ! algal nitrogen
          do iBioTracers = 1, nAlgae
             iIceCount = (iBioTracers-1)*nBioLayersP1
             verticalAlgaeTotalCarbonIceCell(:,iCell) = 0.0_RKIND

             do iLayers = 1,nBioLayersP1
                iBiocount = iBiocount + 1
                verticalAlgaeConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_algaeConc(iBioTracers)+iLayers-1)
                verticalAlgaeIceCell(iLayers+iIceCount,iCell) = verticalAlgaeConcCell(iBioCount,iCell)
                verticalAlgaeTotalCarbonIceCell(iLayers,iCell) = verticalAlgaeTotalCarbonIceCell(iLayers,iCell) + &
                   carbonToNitrogenRatioAlgae(iBioTracers) * verticalAlgaeConcCell(iBioCount,iCell)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBiocount = iBiocount + 1
                verticalAlgaeConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_algaeConc(iBioTracers)+iLayers-1)
             enddo
          enddo
       endif

       ! nitrate
       if (config_use_nitrate) then
          do iLayers = 1, nBioLayersP1
             verticalNitrateConcCell(iLayers,iCell) = tracerArrayCell(tracerObject % index_nitrateConc + iLayers-1)
             verticalNitrateIceCell(iLayers,iCell) = verticalNitrateConcCell(iLayers,iCell)
          enddo
          do iLayers = nBioLayersP1+1, nBioLayersP3
             verticalNitrateConcCell(iLayers,iCell) = tracerArrayCell(tracerObject % index_nitrateConc + iLayers-1)
          enddo
       endif

       if (config_use_carbon) then
          iBioCount = 0
          verticalDOCLabileIceCell(:,iCell) = 0.0_RKIND

          ! DOC
          do iBioTracers = 1, nDOC
             iIceCount = (iBioTracers-1)*nBioLayersP1

             do iLayers = 1,nBioLayersP1
                iBioCount = iBioCount + 1
                verticalDOCConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_DOCConc(iBioTracers) + iLayers-1)
                verticalDOCIceCell(iLayers+iIceCount,iCell) = verticalDOCConcCell(iBioCount,iCell)
                verticalDOCLabileIceCell(iLayers, iCell) = verticalDOCLabileIceCell(iLayers, iCell) + &
                   verticalDOCConcCell(iBioCount,iCell)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDOCConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_DOCConc(iBioTracers) + iLayers-1)
             enddo
          enddo
          iBioCount = 0

          ! DIC
          do iBioTracers = 1, nDIC
             iIceCount = (iBioTracers-1)*nBioLayersP1

             do iLayers = 1,nBioLayersP1
                iBioCount = iBioCount + 1
                verticalDICConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_DICConc(iBioTracers) + iLayers-1)
                verticalDICIceCell(iLayers+iIceCount,iCell) = verticalDICConcCell(iBioCount,iCell)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDICConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_DICConc(iBioTracers) + iLayers-1)
             enddo
          enddo
       endif

       ! DON
       if (config_use_DON) then
          iBioCount = 0
          do iBioTracers = 1, nDON
             iIceCount = (iBioTracers-1)*nBioLayersP1

             do iLayers = 1,nBioLayersP1
                iBioCount = iBioCount + 1
                verticalDONConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_DONConc(iBioTracers) + iLayers-1)
                verticalDONIceCell(iLayers+iIceCount,iCell) = verticalDONConcCell(iBioCount,iCell)
                verticalDOCLabileIceCell(iLayers, iCell) = verticalDOCLabileIceCell(iLayers, iCell) + &
                   verticalDONConcCell(iBioCount,iCell) * carbonToNitrogenRatioDON(iBioTracers)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDONConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_DONConc(iBioTracers) + iLayers-1)
             enddo
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) then
          do iLayers = 1, nBioLayersP1
             verticalAmmoniumConcCell(iLayers,iCell) = &
                  tracerArrayCell(tracerObject % index_ammoniumConc + iLayers-1)
             verticalAmmoniumIceCell(iLayers,iCell) = verticalAmmoniumConcCell(iLayers,iCell)
          enddo
          do iLayers = nBioLayersP1+1, nBioLayersP3
             verticalAmmoniumConcCell(iLayers,iCell) = &
                  tracerArrayCell(tracerObject % index_ammoniumConc + iLayers-1)
          enddo
       endif

       ! silicate
       if (config_use_silicate) then
          do iLayers = 1, nBioLayersP1
             verticalSilicateConcCell(iLayers,iCell) = tracerArrayCell(tracerObject % index_silicateConc+iLayers-1)
             verticalSilicateIceCell(iLayers,iCell) = verticalSilicateConcCell(iLayers,iCell)
          enddo
          do iLayers = nBioLayersP1+1, nBioLayersP3
             verticalSilicateConcCell(iLayers,iCell) = tracerArrayCell(tracerObject % index_silicateConc+iLayers-1)
          enddo
       endif

       ! DMS, DMSPp, DMSPd
       if (config_use_DMS) then
          do iLayers = 1, nBioLayersP1
             verticalDMSConcCell(iLayers,iCell) = tracerArrayCell(tracerObject % index_DMSConc+iLayers-1)
             verticalDMSPpConcCell(iLayers,iCell) = tracerArrayCell(tracerObject % index_DMSPpConc+iLayers-1)
             verticalDMSPdConcCell(iLayers,iCell) = tracerArrayCell(tracerObject % index_DMSPdConc+iLayers-1)
             verticalDMSIceCell(iLayers,iCell)  = verticalDMSConcCell(iLayers,iCell)
             verticalDMSPpIceCell(iLayers,iCell)  = verticalDMSPpConcCell(iLayers,iCell)
             verticalDMSPdIceCell(iLayers,iCell)  = verticalDMSPdConcCell(iLayers,iCell)
          enddo
          do iLayers = nBioLayersP1+1, nBioLayersP3
             verticalDMSConcCell(iLayers,iCell) = tracerArrayCell(tracerObject % index_DMSConc+iLayers-1)
             verticalDMSPpConcCell(iLayers,iCell) = tracerArrayCell(tracerObject % index_DMSPpConc+iLayers-1)
             verticalDMSPdConcCell(iLayers,iCell) = tracerArrayCell(tracerObject % index_DMSPdConc+iLayers-1)
          enddo
       endif

       ! nonreactive tracer
       if (config_use_nonreactive) then
          do iLayers = 1, nBioLayersP1
             verticalNonreactiveConcCell(iLayers,iCell) = &
                  tracerArrayCell(tracerObject % index_nonreactiveConc+iLayers-1)
             verticalNonreactiveIceCell(iLayers,iCell) = verticalNonreactiveConcCell(iLayers,iCell)
          enddo
          do iLayers = nBioLayersP1+1, nBioLayersP3
             verticalNonreactiveConcCell(iLayers,iCell) = &
                  tracerArrayCell(tracerObject % index_nonreactiveConc+iLayers-1)
          enddo
       endif

       ! humic material
       if (config_use_humics) then
          do iLayers = 1, nBioLayersP1
             verticalHumicsConcCell(iLayers,iCell) = tracerArrayCell(tracerObject % index_humicsConc+iLayers-1)
             verticalHumicsIceCell(iLayers,iCell) = verticalHumicsConcCell(iLayers,iCell)
          enddo
          do iLayers = nBioLayersP1+1, nBioLayersP3
             verticalHumicsConcCell(iLayers,iCell) = tracerArrayCell(tracerObject % index_humicsConc+iLayers-1)
          enddo
       endif

       if (config_use_iron) then
          iBioCount = 0

          ! particulate iron
          do iBioTracers = 1, nParticulateIron
             iIceCount = (iBioTracers-1)*nBioLayersP1

             do iLayers = 1,nBioLayersP1
                iBioCount = iBioCount + 1
                verticalParticulateIronConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_particulateIronConc(iBioTracers)+iLayers-1)
                verticalDissolvedIronIceCell(iLayers+iIceCount,iCell) = verticalDissolvedIronConcCell(iBioCount,iCell)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalParticulateIronConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_particulateIronConc(iBioTracers)+iLayers-1)
             enddo
          enddo
          iBioCount = 0

          ! dissolved iron
          do iBioTracers = 1, nDissolvedIron
             iIceCount = (iBioTracers-1)*nBioLayersP1

             do iLayers = 1,nBioLayersP1
                iBioCount = iBioCount + 1
                verticalDissolvedIronConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_dissolvedIronConc(iBioTracers)+iLayers-1)
                verticalDissolvedIronIceCell(iLayers+iIceCount,iCell) = verticalDissolvedIronConcCell(iBioCount,iCell)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDissolvedIronConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_dissolvedIronConc(iBioTracers)+iLayers-1)
             enddo
          enddo
       endif

       ! black carbon and dust aerosols
       if (config_use_zaerosols) then
          iBioCount = 0
          verticalDustTotalIceCell(:,iCell) = 0.0_RKIND
          verticalBCTotalIceCell(:,iCell) = 0.0_RKIND
          verticalDustTotalSnowCell(:,iCell) = 0.0_RKIND
          verticalBCTotalSnowCell(:,iCell) = 0.0_RKIND

          do iBioTracers = 1, nzAerosols
             iIceCount = (iBioTracers-1)*nBioLayersP1
             iSnowCount = (iBioTracers-1)*2
             nAeroType = 1
             if (iBioTracers .gt. 2) nAeroType = 0
             do iLayers = 1,nBioLayersP1
                iBioCount = iBioCount + 1
                verticalAerosolsConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_verticalAerosolsConc(iBioTracers)+iLayers-1)
                verticalAerosolsIceCell(iLayers+iIceCount,iCell) = verticalAerosolsConcCell(iBioCount,iCell)
                verticalBCTotalIceCell(iLayers, iCell) = verticalBCTotalIceCell(iLayers, iCell) + nAeroType * &
                      verticalAerosolsConcCell(iBioCount,iCell)
                verticalDustTotalIceCell(iLayers, iCell) = verticalDustTotalIceCell(iLayers, iCell) + (1-nAeroType) * &
                      verticalAerosolsConcCell(iBioCount,iCell)
             enddo
             do iLayers = nBioLayersP1+1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalAerosolsConcCell(iBioCount,iCell) = &
                     tracerArrayCell(tracerObject % index_verticalAerosolsConc(iBioTracers)+iLayers-1)
                verticalAerosolsSnowCell(iLayers-nBioLayersP1+iSnowCount,iCell) = &
                     verticalAerosolsConcCell(iBioCount,iCell)
                verticalBCTotalSnowCell(iLayers-nBioLayersP1, iCell) = verticalBCTotalSnowCell(iLayers-nBioLayersP1, iCell) + &
                      nAeroType * verticalAerosolsConcCell(iBioCount,iCell)
                verticalDustTotalSnowCell(iLayers-nBioLayersP1,iCell) = verticalDustTotalSnowCell(iLayers-nBioLayersP1,iCell) + &
                      (1-nAeroType) * verticalAerosolsConcCell(iBioCount,iCell)
             enddo
          enddo
       endif

       ! salinity used with BL99 thermodynamics        !echmod deprecate
!       if (config_use_vertical_zsalinity) then
!          do iLayers = 1, nBioLayers
!             verticalSalinityCell(iLayers,iCell) = tracerArrayCell(tracerObject % index_verticalSalinity+iLayers-1)
!          enddo
!       endif
    endif

  end subroutine get_cice_biogeochemistry_tracer_array_cell

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  define_total_biogeochemistry_array_cell
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 16rd September 2024
!> \details
!>
!   Defines total ice and snow arrays for biogeochemistry and aerosols diagnostics
!
!-----------------------------------------------------------------------

  subroutine define_total_biogeochemistry_array_cell(block, tracerObject, totalVerticalBioIceCell, &
         totalVerticalBioSnowCell, iCell)

    type(block_type), intent(in) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    real(kind=RKIND), dimension(:), intent(in):: &
         totalVerticalBioIceCell, &
         totalVerticalBioSnowCell

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_vertical_biochemistry, &
         config_use_nitrate, &
         config_use_carbon, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_nonreactive, &
         config_use_humics, &
         config_use_DON, &
         config_use_iron, &
         config_use_zaerosols

    integer, pointer :: &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON, &
         nParticulateIron, &
         nDissolvedIron, &
         nzAerosols

    type(MPAS_pool_type), pointer :: &
         biogeochemistry

    real(kind=RKIND), dimension(:), pointer :: &
         totalVerticalDiatomIce, &
         totalVerticalSmallPlanktonIce, &
         totalVerticalPhaeocystisIce, &
         totalVerticalNitrateIce, &
         totalVerticalPolysaccharidsIce, &
         totalVerticalLipidsIce, &
         totalVerticalDICIce, &
         totalVerticalAmmoniumIce, &
         totalVerticalSilicateIce, &
         totalVerticalDMSPpIce, &
         totalVerticalDMSPdIce, &
         totalVerticalDMSIce, &
         totalVerticalNonreactiveIce, &
         totalVerticalHumicsIce, &
         totalVerticalProteinsIce, &
         totalVerticalDissolvedIronIce, &
         totalVerticalParticulateIronIce, &
         totalVerticalDustIce, &
         totalVerticalDustSnow, &
         totalVerticalBC1Ice, &
         totalVerticalBC1Snow, &
         totalVerticalBC2Ice, &
         totalVerticalBC2Snow, &
         totalVerticalBCIce, &
         totalVerticalBCSnow, &
         totalVerticalDissolvedIronSnow, &
         totalVerticalAlgaeCarbonIce, &
         totalVerticalDOCLabileIce, &
         carbonToNitrogenRatioAlgae, &
         carbonToNitrogenRatioDON

    call MPAS_pool_get_config(block % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_nitrate", config_use_nitrate)
    call MPAS_pool_get_config(block % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(block % configs, "config_use_ammonium",config_use_ammonium)
    call MPAS_pool_get_config(block % configs, "config_use_silicate",config_use_silicate)
    call MPAS_pool_get_config(block % configs, "config_use_DMS",config_use_DMS)
    call MPAS_pool_get_config(block % configs, "config_use_nonreactive",config_use_nonreactive)
    call MPAS_pool_get_config(block % configs, "config_use_humics",config_use_humics)
    call MPAS_pool_get_config(block % configs, "config_use_DON",config_use_DON)
    call MPAS_pool_get_config(block % configs, "config_use_iron",config_use_iron)
    call MPAS_pool_get_config(block % configs, "config_use_zaerosols",config_use_zaerosols)

    call MPAS_pool_get_dimension(block % dimensions, "nzAerosols", nzAerosols)
    call MPAS_pool_get_dimension(block % dimensions, "nAlgae", nAlgae)
    call MPAS_pool_get_dimension(block % dimensions, "nDOC", nDOC)
    call MPAS_pool_get_dimension(block % dimensions, "nDIC", nDIC)
    call MPAS_pool_get_dimension(block % dimensions, "nDON", nDON)
    call MPAS_pool_get_dimension(block % dimensions, "nParticulateIron", nParticulateIron)
    call MPAS_pool_get_dimension(block % dimensions, "nDissolvedIron", nDissolvedIron)

    call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)

    call MPAS_pool_get_array(biogeochemistry, "totalVerticalDiatomIce", totalVerticalDiatomIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalSmallPlanktonIce", totalVerticalSmallPlanktonIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalPhaeocystisIce", totalVerticalPhaeocystisIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalPolysaccharidsIce", totalVerticalPolysaccharidsIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalLipidsIce", totalVerticalLipidsIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalDICIce", totalVerticalDICIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalProteinsIce", totalVerticalProteinsIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalNitrateIce", totalVerticalNitrateIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalSilicateIce", totalVerticalSilicateIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalAmmoniumIce", totalVerticalAmmoniumIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalDMSIce", totalVerticalDMSIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalDMSPpIce", totalVerticalDMSPpIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalDMSPdIce", totalVerticalDMSPdIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalNonreactiveIce", totalVerticalNonreactiveIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalHumicsIce", totalVerticalHumicsIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalParticulateIronIce", totalVerticalParticulateIronIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalDissolvedIronIce", totalVerticalDissolvedIronIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalDissolvedIronSnow", totalVerticalDissolvedIronSnow)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalDustIce", totalVerticalDustIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalDustSnow", totalVerticalDustSnow)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalBC1Ice", totalVerticalBC1Ice)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalBC1Snow", totalVerticalBC1Snow)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalBC2Ice", totalVerticalBC2Ice)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalBC2Snow", totalVerticalBC2Snow)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalBCIce", totalVerticalBCIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalBCSnow", totalVerticalBCSnow)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalAlgaeCarbonIce", totalVerticalAlgaeCarbonIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalDOCLabileIce", totalVerticalDOCLabileIce)
    call MPAS_pool_get_array(biogeochemistry, "carbonToNitrogenRatioAlgae", carbonToNitrogenRatioAlgae)
    call MPAS_pool_get_array(biogeochemistry, "carbonToNitrogenRatioDON", carbonToNitrogenRatioDON)

    ! biogeochemistry

    if (config_use_vertical_biochemistry) then

       ! algal nitrogen
       totalVerticalDiatomIce(iCell) = totalVerticalBioIceCell(tracerObject % index_algaeConcLayer(1))
       totalVerticalAlgaeCarbonIce(iCell) = totalVerticalDiatomIce(iCell) * carbonToNitrogenRatioAlgae(1)

       SELECT CASE (nAlgae)
       CASE (2)
          totalVerticalSmallPlanktonIce(iCell) = totalVerticalBioIceCell(tracerObject % index_algaeConcLayer(2))
          totalVerticalAlgaeCarbonIce(iCell) = totalVerticalAlgaeCarbonIce(iCell) + &
             totalVerticalSmallPlanktonIce(iCell) * carbonToNitrogenRatioAlgae(2)
       CASE (3)
          totalVerticalSmallPlanktonIce(iCell) = totalVerticalBioIceCell(tracerObject % index_algaeConcLayer(2))
          totalVerticalPhaeocystisIce(iCell) = totalVerticalBioIceCell(tracerObject % index_algaeConcLayer(3))
          totalVerticalAlgaeCarbonIce(iCell) = totalVerticalAlgaeCarbonIce(iCell) + &
             totalVerticalSmallPlanktonIce(iCell) * carbonToNitrogenRatioAlgae(2) + &
             totalVerticalPhaeocystisIce(iCell) * carbonToNitrogenRatioAlgae(3)
       END SELECT
    end if

    ! nitrate
    if (config_use_nitrate) &
      totalVerticalNitrateIce(iCell) = totalVerticalBioIceCell(tracerObject % index_nitrateConcLayer)

    if (config_use_carbon) then

       ! DOC
       totalVerticalDOCLabileIce(iCell) = 0.0_RKIND
       SELECT CASE (nDOC)
       CASE (1)
          totalVerticalPolysaccharidsIce(iCell) = totalVerticalBioIceCell(tracerObject % index_DOCConcLayer(1))
          totalVerticalDOCLabileIce(iCell) = totalVerticalPolysaccharidsIce(iCell)
       CASE (2)
          totalVerticalPolysaccharidsIce(iCell) = totalVerticalBioIceCell(tracerObject % index_DOCConcLayer(1))
          totalVerticalLipidsIce(iCell) = totalVerticalBioIceCell(tracerObject % index_DOCConcLayer(1))
          totalVerticalDOCLabileIce(iCell) = totalVerticalPolysaccharidsIce(iCell) + totalVerticalLipidsIce(iCell)
       END SELECT

       ! DIC
       SELECT CASE (nDIC)
       CASE (1)
          totalVerticalDICIce(iCell) = totalVerticalBioIceCell(tracerObject % index_DICConcLayer(1))
       END SELECT
    end if

    if (config_use_DON) then

       ! DON
       SELECT CASE (nDON)
       CASE (1)
          totalVerticalProteinsIce(iCell) = totalVerticalBioIceCell(tracerObject % index_DONConcLayer(1))
          totalVerticalDOCLabileIce(iCell) = totalVerticalDOCLabileIce(iCell) + &
             totalVerticalProteinsIce(iCell)* carbonToNitrogenRatioDON(1)
       END SELECT
    end if

    ! ammonium
    if (config_use_ammonium) &
      totalVerticalAmmoniumIce(iCell) = totalVerticalBioIceCell(tracerObject % index_ammoniumConcLayer)

    ! silicate
    if (config_use_silicate) &
      totalVerticalSilicateIce(iCell) = totalVerticalBioIceCell(tracerObject % index_silicateConcLayer)

    ! DMS, DMSPp, DMSPd
    if (config_use_DMS) then
       totalVerticalDMSIce(iCell) = totalVerticalBioIceCell(tracerObject % index_DMSConcLayer)
       totalVerticalDMSPpIce(iCell) = totalVerticalBioIceCell(tracerObject % index_DMSPpConcLayer)
       totalVerticalDMSPdIce(iCell) = totalVerticalBioIceCell(tracerObject % index_DMSPdConcLayer)
    end if

    ! nonreactive
    if (config_use_nonreactive) &
      totalVerticalNonreactiveIce(iCell) = totalVerticalBioIceCell(tracerObject % index_nonreactiveConcLayer)

    ! humic material (refractory DOC)
    if (config_use_humics) then
      totalVerticalHumicsIce(iCell) = totalVerticalBioIceCell(tracerObject % index_humicsConcLayer)
    end if

    if (config_use_iron) then

      ! ParticulateIron
      SELECT CASE (nParticulateIron)
      CASE (1)
         totalVerticalParticulateIronIce(iCell) = totalVerticalBioIceCell(tracerObject % index_particulateIronConcLayer(1))
      END SELECT

      ! DissolvedIron
      SELECT CASE (nDissolvedIron)
      CASE (1)
         totalVerticalDissolvedIronIce(iCell) = totalVerticalBioIceCell(tracerObject % index_dissolvedIronConcLayer(1))
         totalVerticalDissolvedIronSnow(iCell) = totalVerticalBioSnowCell(tracerObject % index_dissolvedIronConcLayer(1))
      END SELECT

    end if

    ! black carbon and dust aerosols
    if (config_use_zaerosols) then

      SELECT CASE (nzAerosols)
      CASE (1)
         totalVerticalBC1Ice(iCell) = totalVerticalBioIceCell(tracerObject % index_verticalAerosolsConcLayer(1))
         totalVerticalBC1Snow(iCell) = totalVerticalBioSnowCell(tracerObject % index_verticalAerosolsConcLayer(1))
         totalVerticalBCIce(iCell) = totalVerticalBC1Ice(iCell)
         totalVerticalBCSnow(iCell) = totalVerticalBC1Snow(iCell)
      CASE (2)
         totalVerticalBC1Ice(iCell) = totalVerticalBioIceCell(tracerObject % index_verticalAerosolsConcLayer(1))
         totalVerticalBC1Snow(iCell) = totalVerticalBioSnowCell(tracerObject % index_verticalAerosolsConcLayer(1))
         totalVerticalBC2Ice(iCell) = totalVerticalBioIceCell(tracerObject % index_verticalAerosolsConcLayer(2))
         totalVerticalBC2Snow(iCell) = totalVerticalBioSnowCell(tracerObject % index_verticalAerosolsConcLayer(2))
         totalVerticalBCIce(iCell) = totalVerticalBC1Ice(iCell) + totalVerticalBC2Ice(iCell)
         totalVerticalBCSnow(iCell) = totalVerticalBC1Snow(iCell) + totalVerticalBC2Snow(iCell)
      CASE (3)
         totalVerticalBC1Ice(iCell) = totalVerticalBioIceCell(tracerObject % index_verticalAerosolsConcLayer(1))
         totalVerticalBC1Snow(iCell) = totalVerticalBioSnowCell(tracerObject % index_verticalAerosolsConcLayer(1))
         totalVerticalBC2Ice(iCell) = totalVerticalBioIceCell(tracerObject % index_verticalAerosolsConcLayer(2))
         totalVerticalBC2Snow(iCell) = totalVerticalBioSnowCell(tracerObject % index_verticalAerosolsConcLayer(2))
         totalVerticalDustIce(iCell) = totalVerticalBioIceCell(tracerObject % index_verticalAerosolsConcLayer(3))
         totalVerticalDustSnow(iCell) = totalVerticalBioSnowCell(tracerObject % index_verticalAerosolsConcLayer(3))
         totalVerticalBCIce(iCell) = totalVerticalBC1Ice(iCell) + totalVerticalBC2Ice(iCell)
         totalVerticalBCSnow(iCell) = totalVerticalBC1Snow(iCell) + totalVerticalBC2Snow(iCell)
      END SELECT
   end if

  end subroutine define_total_biogeochemistry_array_cell

!-----------------------------------------------------------------------
! Init CICE parameters
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_icepack_package_parameters
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_icepack_package_parameters(domain, tracerObject)

    type(domain_type), intent(inout) :: domain

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    ! check column configs
    call check_column_package_configs(domain)

    ! set the tracer flags
    call init_icepack_package_tracer_flags(domain)

    ! set the tracer numbers
    call init_icepack_package_tracer_sizes(domain, tracerObject)

    ! set the tracers indices
    call init_icepack_package_tracer_indices(tracerObject)

    ! set the column parameters
    call init_icepack_package_configs(domain)

  end subroutine init_icepack_package_parameters

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  check_column_package_configs
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine check_column_package_configs(domain)

    use seaice_constants, only: &
         seaicePuny

    type(domain_type), intent(inout) :: &
         domain

    integer, pointer :: &
         nCategories, &
         nSnowLayers, &
         nIceLayers, &
         config_nSnowLayers, &
         config_nIceLayers

    character(len=strKIND), pointer :: &
         config_thermodynamics_type, &
         config_heat_conductivity_type, &
         config_shortwave_type, &
         config_albedo_type, &
         config_ice_strength_formulation, &
         config_ridging_participation_function, &
         config_ridging_redistribution_function, &
         config_atmos_boundary_method, &
         config_itd_conversion_type, &
         config_category_bounds_type, &
         config_pond_refreezing_type, &
         config_salt_flux_coupling_type, &
         config_ocean_heat_transfer_type, &
         config_sea_freezing_temperature_type, &
         config_snow_redistribution_scheme

    logical, pointer :: &
         config_calc_surface_stresses, &
         config_calc_surface_temperature, &
         config_use_form_drag, &
         config_use_level_ice, &
         config_use_cesm_meltponds, &  ! deprecated
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_include_pond_freshwater_feedback, &
         config_use_vertical_zsalinity, &      ! deprecated
         config_use_brine, &
         config_use_vertical_tracers, &
         config_use_vertical_biochemistry, &
         config_use_skeletal_biochemistry, &
         config_use_zaerosols, &
         config_use_shortwave_bioabsorption, &
         config_use_nitrate, &
         config_use_carbon, &
         config_use_chlorophyll, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_nonreactive, &
         config_use_humics, &
         config_use_DON, &
         config_use_iron, &
         config_use_modal_aerosols, &
         config_use_column_biogeochemistry, &
         config_use_column_snow_tracers, &
         config_use_effective_snow_density, &
         config_use_snow_liquid_ponds, &
         config_use_snow_grain_radius

    logical :: &
         use_meltponds

    integer :: &
         nPondSchemesActive

    real(kind=RKIND), pointer :: &
         config_max_meltwater_retained_fraction, &
         config_min_meltwater_retained_fraction, &
         config_snow_to_ice_transition_depth

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nCategories", nCategories)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nSnowLayers", nSnowLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nIceLayers", nIceLayers)

    call MPAS_pool_get_config(domain % configs, "config_thermodynamics_type", config_thermodynamics_type)
    call MPAS_pool_get_config(domain % configs, "config_heat_conductivity_type", config_heat_conductivity_type)
    call MPAS_pool_get_config(domain % configs, "config_shortwave_type", config_shortwave_type)
    call MPAS_pool_get_config(domain % configs, "config_albedo_type", config_albedo_type)
    call MPAS_pool_get_config(domain % configs, "config_ice_strength_formulation", config_ice_strength_formulation)
    call MPAS_pool_get_config(domain % configs, "config_ridging_participation_function", config_ridging_participation_function)
    call MPAS_pool_get_config(domain % configs, "config_ridging_redistribution_function", config_ridging_redistribution_function)
    call MPAS_pool_get_config(domain % configs, "config_atmos_boundary_method", config_atmos_boundary_method)
    call MPAS_pool_get_config(domain % configs, "config_itd_conversion_type", config_itd_conversion_type)
    call MPAS_pool_get_config(domain % configs, "config_category_bounds_type", config_category_bounds_type)
    call MPAS_pool_get_config(domain % configs, "config_pond_refreezing_type", config_pond_refreezing_type)
    call MPAS_pool_get_config(domain % configs, "config_salt_flux_coupling_type", config_salt_flux_coupling_type)
    call MPAS_pool_get_config(domain % configs, "config_calc_surface_stresses", config_calc_surface_stresses)
    call MPAS_pool_get_config(domain % configs, "config_calc_surface_temperature", config_calc_surface_temperature)
    call MPAS_pool_get_config(domain % configs, "config_max_meltwater_retained_fraction", config_max_meltwater_retained_fraction)
    call MPAS_pool_get_config(domain % configs, "config_min_meltwater_retained_fraction", config_min_meltwater_retained_fraction)
    call MPAS_pool_get_config(domain % configs, "config_snow_to_ice_transition_depth", config_snow_to_ice_transition_depth)
    call MPAS_pool_get_config(domain % configs, "config_use_form_drag", config_use_form_drag)
    call MPAS_pool_get_config(domain % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds) ! deprecated
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_include_pond_freshwater_feedback", config_include_pond_freshwater_feedback)
    call MPAS_pool_get_config(domain % configs, "config_ocean_heat_transfer_type", config_ocean_heat_transfer_type)
    call MPAS_pool_get_config(domain % configs, "config_sea_freezing_temperature_type", config_sea_freezing_temperature_type)
    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)      ! deprecated
    call MPAS_pool_get_config(domain % configs, "config_use_shortwave_bioabsorption", config_use_shortwave_bioabsorption)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_tracers", config_use_vertical_tracers)
    call MPAS_pool_get_config(domain % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_nitrate", config_use_nitrate)
    call MPAS_pool_get_config(domain % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(domain % configs, "config_use_chlorophyll", config_use_chlorophyll)
    call MPAS_pool_get_config(domain % configs, "config_use_ammonium", config_use_ammonium)
    call MPAS_pool_get_config(domain % configs, "config_use_silicate", config_use_silicate)
    call MPAS_pool_get_config(domain % configs, "config_use_DMS", config_use_DMS)
    call MPAS_pool_get_config(domain % configs, "config_use_nonreactive", config_use_nonreactive)
    call MPAS_pool_get_config(domain % configs, "config_use_humics", config_use_humics)
    call MPAS_pool_get_config(domain % configs, "config_use_DON", config_use_DON)
    call MPAS_pool_get_config(domain % configs, "config_use_iron", config_use_iron)
    call MPAS_pool_get_config(domain % configs, "config_use_modal_aerosols", config_use_modal_aerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_zaerosols", config_use_zaerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
    call MPAS_pool_get_config(domain % configs, "config_nSnowLayers", config_nSnowLayers)
    call MPAS_pool_get_config(domain % configs, "config_nIceLayers", config_nIceLayers)
    call MPAS_pool_get_config(domain % configs, "config_use_column_snow_tracers", config_use_column_snow_tracers)
    call MPAS_pool_get_config(domain % configs, "config_use_effective_snow_density", config_use_effective_snow_density)
    call MPAS_pool_get_config(domain % configs, "config_snow_redistribution_scheme", config_snow_redistribution_scheme)
    call MPAS_pool_get_config(domain % configs, "config_use_snow_liquid_ponds", config_use_snow_liquid_ponds)
    call MPAS_pool_get_config(domain % configs, "config_use_snow_grain_radius", config_use_snow_grain_radius)

    !-----------------------------------------------------------------------
    ! Check values
    !-----------------------------------------------------------------------

    ! deprecate cesm ponds
    if (config_use_cesm_meltponds) then
       call mpas_log_write(&
            "check_column_package_configs: config_use_cesm_meltponds = .true. but cesm ponds have been deprecated", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! deprecate vertical zSalinity
    if (config_use_vertical_zsalinity) then
       call mpas_log_write(&
            "check_column_package_configs: config_use_vertical_zsalinity = .true. but vertical zSalinity has been deprecated", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! deprecate 0-layer thermo
    if (trim(config_thermodynamics_type(1:4)) == "zero") then
       call mpas_log_write(&
            "check_column_package_configs: config_thermodynamics_type) = zero layer but 0-layer thermo has been deprecated", &
            messageType=MPAS_LOG_WARN)
    endif

    ! check config_thermodynamics_type value
    if (.not. (trim(config_thermodynamics_type) == "BL99" .or. &
               trim(config_thermodynamics_type) == "mushy")) then
       call config_error("config_thermodynamics_type", config_thermodynamics_type, "'BL99' or 'mushy'")
    endif

    ! check config_heat_conductivity_type value
    if (.not. (trim(config_heat_conductivity_type) == "MU71" .or. &
               trim(config_heat_conductivity_type) == "bubbly")) then
       call config_error("config_heat_conductivity_type", config_heat_conductivity_type, "'MU71' or 'bubbly'")
    endif

    ! check config_shortwave_type value
    if (.not. (trim(config_shortwave_type) == "ccsm3" .or. &
               trim(config_shortwave_type(1:4)) == "dEdd")) then
       call config_error("config_shortwave_type", config_shortwave_type, "'ccsm3' or 'dEdd'")
    endif

    ! check config_albedo_type value
    if (.not. (trim(config_albedo_type) == "ccsm3" .or. &
               trim(config_albedo_type) == "constant")) then
       call config_error("config_albedo_type", config_albedo_type, "'ccsm3' or 'constant'")
    endif

    ! check config_ice_strength_formulation value
    if (.not. (trim(config_ice_strength_formulation) == "Hibler79" .or. &
               trim(config_ice_strength_formulation) == "Rothrock75")) then
       call config_error("config_ice_strength_formulation", config_ice_strength_formulation, "'Hibler79' or 'Rothrock75'")
    endif

    ! check config_ridging_participation_function value
    if (.not. (trim(config_ridging_participation_function) == "Thorndike75" .or. &
               trim(config_ridging_participation_function) == "exponential")) then
       call config_error("config_ridging_participation_function", &
            config_ridging_participation_function, "'Thorndike75' or 'exponential'")
    endif

    ! check config_ridging_redistribution_function value
    if (.not. (trim(config_ridging_redistribution_function) == "Hibler80" .or. &
               trim(config_ridging_redistribution_function) == "exponential")) then
       call config_error("config_ridging_redistribution_function", &
            config_ridging_redistribution_function, "'Hibler80' or 'exponential'")
    endif

    ! check config_atmos_boundary_method value
    if (.not. (trim(config_atmos_boundary_method) == "similarity" .or. &
               trim(config_atmos_boundary_method) == "constant" .or. &
               trim(config_atmos_boundary_method) == "mixed")) then
       call config_error("config_atmos_boundary_method", config_atmos_boundary_method, "'similarity' or 'constant' or 'mixed'")
    endif

    ! check config_itd_conversion_type value
    if (.not. (trim(config_itd_conversion_type) == "delta function" .or. &
               trim(config_itd_conversion_type) == "linear remap")) then
       call config_error("config_itd_conversion_type", config_itd_conversion_type, "'delta function' or 'linear remap'")
    endif

    ! check config_category_bounds_type value
    if (.not. (trim(config_category_bounds_type) == "single category" .or. &
               trim(config_category_bounds_type) == "original" .or. &
               trim(config_category_bounds_type) == "new" .or. &
               trim(config_category_bounds_type) == "WMO" .or. &
               trim(config_category_bounds_type) == "asymptotic")) then
       call config_error("config_category_bounds_type", &
            config_category_bounds_type, "'single category', 'original', 'new', 'WMO' or 'asymptotic'")
    endif

    ! check config_pond_refreezing_type value
    if (.not. (trim(config_pond_refreezing_type) == "cesm" .or. &
               trim(config_pond_refreezing_type) == "hlid")) then
       call config_error("config_pond_refreezing_type", config_pond_refreezing_type, "'cesm' or 'hlid'")
    endif

    ! check config_salt_flux_coupling_type
    if (trim(config_salt_flux_coupling_type) /= "constant") then
       call config_error("config_salt_flux_coupling_type", config_salt_flux_coupling_type, "'constant'")
    endif

    ! check for consistency in snow vertical dimension
    if (config_nSnowLayers /= nSnowlayers) &
       call mpas_log_write(&
            'Check for inconsistencies in restart file: config_nSnowLayers /= nSnowLayers', &
            messageType=MPAS_LOG_CRIT)

    ! check for consistency in ice vertical dimension
    if (config_nIceLayers /= nIcelayers) &
       call mpas_log_write(&
            'Check for inconsistencies in restart file: config_nIceLayers /= nIceLayers', &
            messageType=MPAS_LOG_CRIT)

    !-----------------------------------------------------------------------
    ! Check combinations
    !-----------------------------------------------------------------------

    ! check only single meltpond option on
    nPondSchemesActive = 0
    if (config_use_level_meltponds) nPondSchemesActive = nPondSchemesActive + 1
    if (config_use_topo_meltponds)  nPondSchemesActive = nPondSchemesActive + 1
    if (nPondSchemesActive > 1) then
       call mpas_log_write(&
            'check_column_package_configs: More than one melt pond scheme active', &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check for itd remapping with only one category
    if (nCategories == 1 .and. trim(config_itd_conversion_type) == "linear remap") then
       call mpas_log_write(&
            'check_column_package_configs: Remapping the ITD is not allowed for nCategories=1', &
            messageType=MPAS_LOG_ERR)
       call mpas_log_write(&
            "Use config_itd_conversion_type = 'delta function' with config_category_bounds_type = 'original'", &
            messageType=MPAS_LOG_ERR)
       call mpas_log_write(&
            "or for column configurations use config_category_bounds_type = 'single category'", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check itd and category bounds discrepancy
    if (nCategories /= 1 .and. trim(config_category_bounds_type) == 'single category') then
       call mpas_log_write(&
            "check_column_package_configs: nCategories /= 1 .and. config_category_bounds_type = 'single category'", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check level ponds and level ice
    if (config_use_level_meltponds .and. .not. config_use_level_ice) then
       call mpas_log_write(&
            "check_column_package_configs: config_use_level_meltponds = .true. and config_use_level_ice = .false.", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check config_snow_to_ice_transition_depth and level ponds
    if (config_use_level_meltponds .and. abs(config_snow_to_ice_transition_depth) > seaicePuny) then
       call mpas_log_write(&
            "check_column_package_configs: config_use_level_meltponds = .true. and config_snow_to_ice_transition_depth /= 0", &
            messageType=MPAS_LOG_CRIT)
    endif

    if (config_use_level_meltponds .and. config_include_pond_freshwater_feedback) then
       call mpas_log_write(&
            "check_column_package_configs: config_include_pond_freshwater_feedback is only available for topo ponds", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check config_snow_redistribution_scheme
    if (trim(config_snow_redistribution_scheme) == "ITDrdg".and. .not. config_use_effective_snow_density) then
       call mpas_log_write(&
            "check_column_package_configs: config_snow_redistribution = 'ITDrdg' but config_use_effective_snow_density is false", &
            messageType=MPAS_LOG_CRIT)
    endif
    if (trim(config_snow_redistribution_scheme) == "ITDrdg".and. .not. config_use_level_ice) then
       call mpas_log_write(&
            "check_column_package_configs: config_snow_redistribution = 'ITDrdg' but config_use_level_ice = .false.", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check config_use_snow_liquid_ponds and config_use_snow_grain_radius
    if (config_use_snow_liquid_ponds .and. .not. config_use_snow_grain_radius) then
       call mpas_log_write(&
            "check_column_package_configs: config_use_snow_liquid_ponds = true but config_use_snow_grain_radius = false", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check snow grain radius and dEdd shortwave
    if (config_use_snow_grain_radius .and. .not. (trim(config_shortwave_type(1:4)) == "dEdd")) then
       call mpas_log_write(&
            "check_column_package_configs: snow grain radius requires config_shortwave_type(1:4) ==  'dEdd'", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check dEdd shortwave if using ponds
    use_meltponds = (config_use_level_meltponds .or. config_use_topo_meltponds)
    if (trim(config_shortwave_type(1:4)) /= 'dEdd' .and. use_meltponds .and. config_calc_surface_temperature) then
       call mpas_log_write(&
            "check_column_package_configs: config_shortwave_type(1:4)) /= 'dEdd' .and. use_meltponds = .true.", &
            messageType=MPAS_LOG_ERR)
       call mpas_log_write(&
            ".and. config_calc_surface_temperature ==.true.", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check range of config_min_meltwater_retained_fraction and config_max_meltwater_retained_fraction
    if (config_min_meltwater_retained_fraction < 0.0_RKIND .or. &
        config_min_meltwater_retained_fraction > 1.0_RKIND) then
       call mpas_log_write(&
            'check_column_package_configs: config_min_meltwater_retained_fraction out of bounds', &
            messageType=MPAS_LOG_CRIT)
    endif
    if (config_max_meltwater_retained_fraction < 0.0_RKIND .or. &
        config_max_meltwater_retained_fraction > 1.0_RKIND) then
       call mpas_log_write(&
            'check_column_package_configs: config_max_meltwater_retained_fraction out of bounds', &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check mushy physics and surface temperature calculation
    if (trim(config_thermodynamics_type) == "mushy" .and. .not. config_calc_surface_temperature) then
       call mpas_log_write(&
            "check_column_package_configs: config_thermodynamics_type = 'mushy' and config_calc_surface_temperature = .false.", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check not form drag with constant atmosphere boundary method
    if (config_use_form_drag .and. trim(config_atmos_boundary_method) == "constant") then
       call mpas_log_write(&
            "check_column_package_configs: config_use_form_drag = .true. and config_atmos_boundary_method = 'constant'", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check not form drag with not calculating surface stresses
    if (config_use_form_drag .and. .not. config_calc_surface_stresses) then
       call mpas_log_write(&
            "check_column_package_configs: config_use_form_drag = .true. and config_calc_surface_stresses = .false.", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check using form drag but not level ice
    if (config_use_form_drag .and. .not. config_use_level_ice) then
       call mpas_log_write(&
            "check_column_package_configs: config_use_form_drag = .true. and config_use_level_ice = .false.", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check form drag and ocean heat flux type
    if (.not. config_use_form_drag .and. trim(config_ocean_heat_transfer_type) == "Cdn_ocn") then
       call mpas_log_write(&
            "check_column_package_configs: config_use_form_drag = .false. and config_ocean_heat_transfer_type == 'Cdn_ocn'", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check thermodynamic type and sea freezing temperature type
    if (trim(config_thermodynamics_type) == "BL99" .and. trim(config_sea_freezing_temperature_type) /= "linear_salt") then
       call mpas_log_write(&
            "check_column_package_configs: config_thermodynamics_type == 'BL99' "//&
            "and config_sea_freezing_temperature_type /= 'linear_salt'", &
            messageType=MPAS_LOG_CRIT)
    endif
    if (trim(config_thermodynamics_type) == "mushy" .and. trim(config_sea_freezing_temperature_type) /= "mushy") then
       call mpas_log_write(&
            "check_column_package_configs: config_thermodynamics_type == 'mushy' and "//&
            "config_sea_freezing_temperature_type /= 'mushy'", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check biogeochemistry flags:

    if (.not. config_use_column_biogeochemistry .and. (config_use_vertical_biochemistry .or. &
        config_use_skeletal_biochemistry .or. config_use_nitrate .or. config_use_carbon .or. &
        config_use_chlorophyll .or. config_use_ammonium .or. config_use_silicate .or. &
        config_use_DMS .or. config_use_nonreactive .or. config_use_humics .or. &
        config_use_DON .or. config_use_iron)) then

       call mpas_log_write(&
            "check_column_package_configs: config_use_column_biogeochemistry = false. "//&
            "All biogeochemistry namelist flags must also be false", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check that vertical bio tracers use brine height
    if ((config_use_vertical_biochemistry .or. config_use_zaerosols) .and. &
         (.not. config_use_brine .or. .not. config_use_vertical_tracers )) then
       call mpas_log_write(&
            "check_column_package_configs: vertical biochemistry and zaerosols require " //&
            "config_use_brine and config_use_vertical_tracer = true", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check that vertical bio tracers and not used with skeletal bio tracers
    if (config_use_vertical_tracers  .and.  config_use_skeletal_biochemistry) then
       call mpas_log_write(&
            "check_column_package_configs: vertical bio tracers and skeletal bio tracers cannot both be true", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check that brine height is used with either aerosols or bgc
    if (config_use_brine  .and. &
         (.not. config_use_column_biogeochemistry .and. .not. config_use_zaerosols)) then
       call mpas_log_write(&
            "check_column_package_configs: brine tracer must be used with vertical tracers - config_use_column_biogeochemistry and/or config_use_zaerosols equal to true", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check that the shortwave scheme and bioabsorption is consistent
    if (config_use_shortwave_bioabsorption .and. .not. (trim(config_shortwave_type(1:4)) == "dEdd")) then
       call mpas_log_write(&
            "check_column_package_configs: shortwave bioabsorption requires config_shortwave_type(1:4) ==  'dEdd'", &
            messageType=MPAS_LOG_CRIT)
    endif

    ! check that nitrate is true for biogeochemistry
    if ((config_use_vertical_biochemistry .or. config_use_skeletal_biochemistry) .and. .not. config_use_nitrate) then
       call mpas_log_write(&
            "check_column_package_configs: biochemistry needs at the very least config_use_nitrate = true", &
            messageType=MPAS_LOG_CRIT)
    endif

  end subroutine check_column_package_configs

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_icepack_package_tracer_flags
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_icepack_package_tracer_flags(domain)

    use icepack_intfc, only: &
         icepack_init_tracer_flags

    type(domain_type), intent(inout) :: domain

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols, &
         config_use_brine, &
!         config_use_vertical_zsalinity, &      !echmod deprecate
         config_use_zaerosols, &
         config_use_nitrate, &
         config_use_DON, &
         config_use_carbon, &
         config_use_chlorophyll, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_iron, &
         config_use_humics, &
         config_use_nonreactive, &
         config_use_vertical_biochemistry, &
         config_use_skeletal_biochemistry, &
         config_use_effective_snow_density, &
         config_use_snow_grain_radius

    logical :: &
         use_meltponds, &
         use_nitrogen

    call MPAS_pool_get_config(domain % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(domain % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_aerosols", config_use_aerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)
!    call MPAS_pool_get_config(domain % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)   !echmod deprecate
    call MPAS_pool_get_config(domain % configs, "config_use_zaerosols", config_use_zaerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_nitrate", config_use_nitrate)
    call MPAS_pool_get_config(domain % configs, "config_use_DON", config_use_DON)
    call MPAS_pool_get_config(domain % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(domain % configs, "config_use_chlorophyll", config_use_chlorophyll)
    call MPAS_pool_get_config(domain % configs, "config_use_ammonium", config_use_ammonium)
    call MPAS_pool_get_config(domain % configs, "config_use_silicate", config_use_silicate)
    call MPAS_pool_get_config(domain % configs, "config_use_DMS", config_use_DMS)
    call MPAS_pool_get_config(domain % configs, "config_use_iron", config_use_iron)
    call MPAS_pool_get_config(domain % configs, "config_use_humics", config_use_humics)
    call MPAS_pool_get_config(domain % configs, "config_use_nonreactive", config_use_nonreactive)
    call MPAS_pool_get_config(domain % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_effective_snow_density", config_use_effective_snow_density)
    call MPAS_pool_get_config(domain % configs, "config_use_snow_grain_radius", config_use_snow_grain_radius)

    use_nitrogen = .false.
    if (config_use_skeletal_biochemistry .or. config_use_vertical_biochemistry) &
         use_nitrogen = .true.

    use_meltponds = (config_use_level_meltponds .or. config_use_topo_meltponds)

    call icepack_init_tracer_flags(&
         tr_iage_in      = config_use_ice_age, &
         tr_FY_in        = config_use_first_year_ice, &
         tr_lvl_in       = config_use_level_ice, &
         tr_pond_in      = use_meltponds, &
         tr_pond_lvl_in  = config_use_level_meltponds, &
         tr_pond_topo_in = config_use_topo_meltponds, &
         !tr_fsd_in       = config_use_floe_size_distribution, &
         tr_aero_in      = config_use_aerosols, &
         !tr_iso_in       = , &
         tr_brine_in     = config_use_brine, &
         tr_zaero_in     = config_use_zaerosols, &
         tr_bgc_Nit_in   = config_use_nitrate, &
         tr_bgc_N_in     = use_nitrogen, &
         tr_bgc_DON_in   = config_use_DON, &
         tr_bgc_C_in     = config_use_carbon, &
         tr_bgc_chl_in   = config_use_chlorophyll, &
         tr_bgc_Am_in    = config_use_ammonium, &
         tr_bgc_Sil_in   = config_use_silicate, &
         tr_bgc_DMS_in   = config_use_DMS, &
         tr_bgc_Fe_in    = config_use_iron, &
         tr_bgc_hum_in   = config_use_humics, &
         tr_bgc_PON_in   = config_use_nonreactive)

    call seaice_icepack_write_warnings(icepack_warnings_aborted())

  end subroutine init_icepack_package_tracer_flags

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_icepack_package_tracer_sizes
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_icepack_package_tracer_sizes(domain, tracerObject)

    use icepack_intfc, only: &
         icepack_init_tracer_sizes

    type(domain_type), intent(in) :: &
         domain

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    integer, pointer :: &
         nCategories, &
         nIceLayers, &
         nSnowLayers, &
         nAerosols, &
         nBioLayers, &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON, &
         nParticulateIron, &
         nDissolvedIron, &
         nzAerosols

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nCategories", nCategories)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nSnowLayers", nSnowLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nBioLayers", nBioLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nAerosols", nAerosols)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nzAerosols", nzAerosols)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nAlgae", nAlgae)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDOC", nDOC)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDIC", nDIC)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDON", nDON)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nParticulateIron", nParticulateIron)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDissolvedIron", nDissolvedIron)

    call icepack_init_tracer_sizes(&
         ncat_in      = nCategories, &
         nilyr_in     = nIceLayers, &
         nslyr_in     = nSnowLayers, &
         nblyr_in     = nBioLayers, &
         !nfsd_in      = , &
         !n_iso_in     = , &
         !n_aero_in    = nAerosols, &
         n_algae_in   = nAlgae, &
         n_DOC_in     = nDOC, &
         n_DON_in     = nDON, &
         n_DIC_in     = nDIC, &
         n_fed_in     = nDissolvedIron, &
         n_fep_in     = nParticulateIron, &
         n_zaero_in   = nzAerosols, &
         ntrcr_in     = tracerObject % nTracers, &
         ntrcr_o_in   = tracerObject % nTracersNotBio, &
         nbtrcr_in    = tracerObject % nBioTracers, &
         nbtrcr_sw_in = tracerObject % nBioTracersShortwave)

    call seaice_icepack_write_warnings(icepack_warnings_aborted())

  end subroutine init_icepack_package_tracer_sizes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_icepack_package_tracer_indices
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_icepack_package_tracer_indices(tracerObject)

    use icepack_intfc, only: &
         icepack_init_tracer_indices

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    call icepack_init_tracer_indices(&
         nt_Tsfc_in       = tracerObject % index_surfaceTemperature, &
         nt_qice_in       = tracerObject % index_iceEnthalpy, &
         nt_qsno_in       = tracerObject % index_snowEnthalpy, &
         nt_sice_in       = tracerObject % index_iceSalinity, &
         nt_fbri_in       = tracerObject % index_brineFraction, &
         nt_iage_in       = tracerObject % index_iceAge, &
         nt_FY_in         = tracerObject % index_firstYearIceArea, &
         nt_alvl_in       = tracerObject % index_levelIceArea, &
         nt_vlvl_in       = tracerObject % index_levelIceVolume, &
         nt_apnd_in       = tracerObject % index_pondArea, &
         nt_hpnd_in       = tracerObject % index_pondDepth, &
         nt_ipnd_in       = tracerObject % index_pondLidThickness, &
         nt_smice_in      = tracerObject % index_snowIceMass, &
         nt_smliq_in      = tracerObject % index_snowLiquidMass, &
         nt_rhos_in       = tracerObject % index_snowDensity, &
         nt_rsnw_in       = tracerObject % index_snowGrainRadius, &
         !nt_fsd
         !nt_isosno_in
         !nt_isoice_in
         nt_aero_in       = tracerObject % index_aerosols, &
         nt_zaero_in      = tracerObject % index_verticalAerosolsConc, & ! name change?nt_zaeros -> nt_zaero
         !nt_bgc_C_in      = tracerObject % index_algalCarbon, & ! not supported
         nt_bgc_N_in      = tracerObject % index_algaeConc, &
         nt_bgc_chl_in    = tracerObject % index_algalChlorophyll, &
         nt_bgc_DOC_in    = tracerObject % index_DOCConc, &
         nt_bgc_DON_in    = tracerObject % index_DONConc, &
         nt_bgc_DIC_in    = tracerObject % index_DICConc, &
         nt_bgc_Fed_in    = tracerObject % index_dissolvedIronConc, &
         nt_bgc_Fep_in    = tracerObject % index_particulateIronConc, &
         nt_bgc_Nit_in    = tracerObject % index_nitrateConc, &
         nt_bgc_Am_in     = tracerObject % index_ammoniumConc, &
         nt_bgc_Sil_in    = tracerObject % index_silicateConc, &
         nt_bgc_DMSPp_in  = tracerObject % index_DMSPpConc, &
         nt_bgc_DMSPd_in  = tracerObject % index_DMSPdConc, &
         nt_bgc_DMS_in    = tracerObject % index_DMSConc, &
         nt_bgc_hum_in    = tracerObject % index_humicsConc, &
         nt_bgc_PON_in    = tracerObject % index_nonreactiveConc, &
         nlt_zaero_in     = tracerObject % index_verticalAerosolsConcLayer, &
!         nlt_bgc_C_in     = tracerObject % index_algalCarbonLayer, &                !echmod: not fully supported
         nlt_bgc_N_in     = tracerObject % index_algaeConcLayer, &
         nlt_bgc_chl_in   = tracerObject % index_algalChlorophyllLayer, &
         nlt_bgc_DOC_in   = tracerObject % index_DOCConcLayer, &
         nlt_bgc_DON_in   = tracerObject % index_DONConcLayer, &
         nlt_bgc_DIC_in   = tracerObject % index_DICConcLayer, &
         nlt_bgc_Fed_in   = tracerObject % index_dissolvedIronConcLayer, &
         nlt_bgc_Fep_in   = tracerObject % index_particulateIronConcLayer, &
         nlt_bgc_Nit_in   = tracerObject % index_nitrateConcLayer, &
         nlt_bgc_Am_in    = tracerObject % index_ammoniumConcLayer, &
         nlt_bgc_Sil_in   = tracerObject % index_silicateConcLayer, &
         nlt_bgc_DMSPp_in = tracerObject % index_DMSPpConcLayer, &
         nlt_bgc_DMSPd_in = tracerObject % index_DMSPdConcLayer, &
         nlt_bgc_DMS_in   = tracerObject % index_DMSConcLayer, &
         nlt_bgc_hum_in   = tracerObject % index_humicsConcLayer, &
         nlt_bgc_PON_in   = tracerObject % index_nonreactiveConcLayer, &
         nt_zbgc_frac_in  = tracerObject % index_mobileFraction, &
         nt_bgc_S_in      = tracerObject % index_verticalSalinity, & ! name change nt_zbgc_S->nt_bgc_S_in
         nlt_chl_sw_in    = tracerObject % index_chlorophyllShortwave, &
         nlt_zaero_sw_in  = tracerObject % index_verticalAerosolsConcShortwave, &
         bio_index_o_in   = tracerObject % index_LayerIndexToDataArray, &
         bio_index_in     = tracerObject % index_LayerIndexToBioIndex)

    call seaice_icepack_write_warnings(icepack_warnings_aborted())

  end subroutine init_icepack_package_tracer_indices

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_icepack_package_configs
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_icepack_package_configs(domain)

    use icepack_intfc, only: &
         icepack_init_parameters, &
         icepack_write_parameters, &
         icepack_configure, &
         icepack_query_parameters ! debugging

    use seaice_constants, only: &
         seaicePuny, &                       ! a small number
         seaiceBigNumber, &                     ! a large number
         seaiceSecondsPerDay, &              ! number of seconds in 1 day
         seaiceDensityIce, &                 ! density of ice (kg/m^3)
         seaiceDensitySnow, &                ! density of snow (kg/m^3)
         seaiceDensitySeaWater, &            ! density of seawater (kg/m^3)
         seaiceDensityFreshwater, &          ! density of freshwater (kg/m^3)
         seaiceFreshIceSpecificHeat, &       ! specific heat of fresh ice (J/kg/K)
         seaiceWaterVaporSpecificHeat, &     ! specific heat of water vapor (J/kg/K)
         seaiceAirSpecificHeat, &            ! specific heat of air (J/kg/K)
         seaiceSeaWaterSpecificHeat, &       ! specific heat of ocn (J/kg/K)
         seaiceLatentHeatVaporization, &     ! latent heat, vaporization freshwater (J/kg)
         seaiceZvir, &                       ! rh2o/rair - 1.0
         seaiceLatentHeatSublimation, &      ! latent heat, sublimation freshwater (J/kg)
         seaiceLatentHeatMelting, &          ! latent heat of melting of fresh ice (J/kg)
         seaiceIceSurfaceMeltingTemperature, &  ! melting temp. ice top surface  (C)
         seaiceSnowSurfaceMeltingTemperature, & ! melting temp. snow top surface  (C)
         seaiceStefanBoltzmann, &            ! J m-2 K-4 s-1
         seaiceIceSnowEmissivity, &          ! emissivity of snow and ice
         seaiceSnowSurfaceScatteringLayer, & ! snow surface scattering layer thickness (m)
         seaiceStabilityReferenceHeight, &   ! stability reference height (m)
         seaiceFreshWaterFreezingPoint, &    ! freezing temp of fresh ice (K)
         seaiceIceSurfaceRoughness, &        ! ice surface roughness (m)
         seaiceVonKarmanConstant, &          ! Von Karman constant
         seaiceOceanAlbedo, &                ! Ocean albedo
         seaiceReferenceSalinity, &          ! ice reference salinity (ppt)
         seaiceMaximumSalinity, &            ! ice maximum salinity (ppt)
         seaiceMeltingTemperatureDepression, &  ! melting temperature depression factor (C/ppt)
         seaiceFrazilSalinityReduction, &    ! bulk salinity reduction of newly formed frazil (ppt)
         seaiceFrazilIcePorosity, &          ! initial liquid fraction of frazil
         seaicePi, &                         ! pi
         seaiceGravity, &                    ! gravitational acceleration (m/s^2)
         seaiceSnowPatchiness, &             ! snow patchiness parameter
         seaiceIceStrengthConstantHiblerP, & ! P* constant in Hibler strength formulation
         seaiceIceStrengthConstantHiblerC, & ! C* constant in Hibler strength formulation
         skeletalLayerThickness, &           ! skeletal layer thickness
         gramsCarbonPerMolCarbon, &          ! grams carbon per mol carbon
         seaiceSnowMinimumDensity, &         ! minimum snow density (kg/m^3)
         seaiceBrineDynamicViscosity, &      ! dynamic viscosity of brine (kg/m/s)
         seaiceFreezingTemperatureConstant, &! constant freezing temp of seawater (C)
         seaiceExtinctionCoef, &             ! vis extnctn coef in ice, wvlngth<700nm (1/m)
         seaiceFreshIceConductivity, &       ! thermal conductivity of fresh ice(W/m/deg)
         seaiceSnowMinimumThickness, &       ! min snow thickness for computing zTsn (m)
         seaiceFrazilMinimumThickness, &     ! min thickness of new frazil ice
         seaiceAlbedoWtVisibleDirect, &      ! visible, direct
         seaiceAlbedoWtNearIRDirect, &       ! near IR, direct
         seaiceAlbedoWtVisibleDiffuse, &     ! visible, diffuse
         seaiceAlbedoWtNearIRDiffuse, &      ! near IR, diffuse
         seaiceQsatQiceConstant, &           ! constant for saturation humidity over ice
         seaiceQsatTiceConstant, &           ! constant for saturation humidity over ice
         seaiceQsatQocnConstant, &           ! constant for saturation humidity over ocean
         seaiceQsatTocnConstant              ! constant for saturation humidity over ocean

    type(domain_type), intent(inout) :: &
         domain

    type(MPAS_pool_type), pointer :: &
         snow

    character(len=strKIND), pointer :: &
         config_thermodynamics_type, &
         config_heat_conductivity_type, &
         config_shortwave_type, &
         config_albedo_type, &
         config_ice_strength_formulation, &
         config_ridging_participation_function, &
         config_ridging_redistribution_function, &
         config_atmos_boundary_method, &
         config_itd_conversion_type, &
         config_category_bounds_type, &
         config_pond_refreezing_type, &
         config_salt_flux_coupling_type, &
         config_ocean_heat_transfer_type, &
         config_frazil_coupling_type, &
         config_congelation_freezing_method, &
         config_sea_freezing_temperature_type, &
         config_skeletal_bgc_flux_type, &
         config_snow_redistribution_scheme

    logical, pointer :: &
         config_use_snicar_ad, &
         config_calc_surface_temperature, &
         config_update_ocean_fluxes, &
         config_use_form_drag, &
         config_use_high_frequency_coupling, &
         config_use_ocean_mixed_layer, &
         config_calc_surface_stresses, &
         config_use_vertical_tracers, &
         config_scale_initial_vertical_bgc, &
         config_use_vertical_biochemistry, &
         config_use_shortwave_bioabsorption, &
         config_use_skeletal_biochemistry, &
         config_use_modal_aerosols, &
         config_use_macromolecules, &
         config_do_restart_bgc, &
         config_use_snow_liquid_ponds, &
         config_use_snow_grain_radius, &
         config_use_shortwave_redistribution, &
         config_use_iron_solubility_file

    real(kind=RKIND), pointer :: &
         config_min_friction_velocity, &
         config_ice_ocean_drag_coefficient, &
         config_snow_thermal_conductivity, &
         config_rapid_mode_channel_radius, &
         config_rapid_model_critical_Ra, &
         config_rapid_mode_aspect_ratio, &
         config_slow_mode_drainage_strength, &
         config_slow_mode_critical_porosity, &
         config_macro_drainage_timescale, &
         config_congelation_ice_porosity, &
!         config_frazil_ice_porosity, &
!         config_frazil_salinity_reduction, &
         config_visible_ice_albedo, &
         config_infrared_ice_albedo, &
         config_visible_snow_albedo, &
         config_infrared_snow_albedo, &
         config_variable_albedo_thickness_limit, &
!         config_snow_surface_scattering_layer_depth, &
         config_ice_shortwave_tuning_parameter, &
         config_pond_shortwave_tuning_parameter, &
         config_snow_shortwave_tuning_parameter, &
         config_shortwave_redistribution_fraction, &
         config_shortwave_redistribution_threshold, &
         config_temp_change_snow_grain_radius_change, &
         config_max_melting_snow_grain_radius, &
         config_algae_absorption_coefficient, &
         config_ridging_efolding_scale, &
         config_ratio_ridging_work_to_PE, &
         config_snow_to_ice_transition_depth, &
         config_pond_flushing_factor, &
         config_min_meltwater_retained_fraction, &
         config_max_meltwater_retained_fraction, &
         config_pond_depth_to_fraction_ratio, &
         config_snow_on_pond_ice_tapering_parameter, &
         config_critical_pond_ice_thickness, &
         config_biogrid_bottom_molecular_sublayer, &
         config_bio_gravity_drainage_length_scale, &
         config_biogrid_top_molecular_sublayer, &
         config_snow_porosity_at_ice_surface, &
         config_ratio_Si_to_N_diatoms, &
         config_ratio_Si_to_N_small_plankton, &
         config_ratio_Si_to_N_phaeocystis, &
         config_ratio_S_to_N_diatoms, &
         config_ratio_S_to_N_small_plankton, &
         config_ratio_S_to_N_phaeocystis, &
         config_ratio_Fe_to_C_diatoms, &
         config_ratio_Fe_to_C_small_plankton, &
         config_ratio_Fe_to_C_phaeocystis, &
         config_ratio_Fe_to_N_diatoms, &
         config_ratio_Fe_to_N_small_plankton, &
         config_ratio_Fe_to_N_phaeocystis, &
         config_ratio_Fe_to_DON, &
         config_ratio_Fe_to_DOC_saccharids, &
         config_ratio_Fe_to_DOC_lipids, &
         config_respiration_fraction_of_growth, &
         config_algal_maximum_velocity, &
         config_ratio_Fe_to_dust, &
         config_solubility_of_Fe_in_dust, &
         config_chla_absorptivity_of_diatoms, &
         config_chla_absorptivity_of_small_plankton, &
         config_chla_absorptivity_of_phaeocystis, &
         config_light_attenuation_diatoms, &
         config_light_attenuation_small_plankton, &
         config_light_attenuation_phaeocystis, &
         config_light_inhibition_diatoms, &
         config_light_inhibition_small_plankton, &
         config_light_inhibition_phaeocystis, &
         config_maximum_growth_rate_diatoms, &
         config_maximum_growth_rate_small_plankton, &
         config_maximum_growth_rate_phaeocystis, &
         config_temperature_growth_diatoms, &
         config_temperature_growth_small_plankton, &
         config_temperature_growth_phaeocystis, &
         config_grazed_fraction_diatoms, &
         config_grazed_fraction_small_plankton, &
         config_grazed_fraction_phaeocystis, &
         config_mortality_diatoms, &
         config_mortality_small_plankton, &
         config_mortality_phaeocystis, &
         config_temperature_mortality_diatoms, &
         config_temperature_mortality_small_plankton, &
         config_temperature_mortality_phaeocystis, &
         config_exudation_diatoms, &
         config_exudation_small_plankton, &
         config_exudation_phaeocystis, &
         config_nitrate_saturation_diatoms, &
         config_nitrate_saturation_small_plankton, &
         config_nitrate_saturation_phaeocystis, &
         config_ammonium_saturation_diatoms, &
         config_ammonium_saturation_small_plankton, &
         config_ammonium_saturation_phaeocystis, &
         config_silicate_saturation_diatoms, &
         config_silicate_saturation_small_plankton, &
         config_silicate_saturation_phaeocystis, &
         config_iron_saturation_diatoms, &
         config_iron_saturation_small_plankton, &
         config_iron_saturation_phaeocystis, &
         config_fraction_spilled_to_DON, &
         config_degredation_of_DON, &
         config_fraction_DON_ammonium, &
         config_fraction_loss_to_saccharids, &
         config_fraction_loss_to_lipids, &
         config_fraction_exudation_to_saccharids, &
         config_fraction_exudation_to_lipids, &
         config_remineralization_saccharids, &
         config_remineralization_lipids, &
         config_maximum_brine_temperature, &
         config_salinity_dependence_of_growth, &
         config_minimum_optical_depth, &
         config_slopped_grazing_fraction, &
         config_excreted_fraction, &
         config_fraction_mortality_to_ammonium, &
         config_fraction_iron_remineralized, &
         config_nitrification_rate, &
         config_desorption_loss_particulate_iron, &
         config_maximum_loss_fraction, &
         config_maximum_ratio_iron_to_saccharids, &
         config_respiration_loss_to_DMSPd, &
         config_DMSP_to_DMS_conversion_fraction, &
         config_DMSP_to_DMS_conversion_time, &
         config_DMS_oxidation_time, &
         config_fallen_snow_radius, &
         config_new_snow_density, &
         config_max_snow_density, &
         config_minimum_wind_compaction, &
         config_wind_compaction_factor, &
         config_snow_redistribution_factor, &
         config_max_dry_snow_radius, &
         config_floediam, &
         config_floeshape

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         snowEmpiricalGrowthParameterTau, &
         snowEmpiricalGrowthParameterKappa, &
         snowPropertyRate

    integer, pointer :: &
         config_boundary_layer_iteration_number, &
         nGrainAgingTemperature, &
         nGrainAgingTempGradient, &
         nGrainAgingSnowDensity

    integer :: &
         config_thermodynamics_type_int, &
         config_ice_strength_formulation_int, &
         config_ridging_participation_function_int, &
         config_ridging_redistribution_function_int, &
         config_itd_conversion_type_int, &
         config_category_bounds_type_int

    character(len=strKIND) :: tmp_config_shortwave_type
    character(len=strKIND) :: tmp_config_snw_ssp_table
    character(len=strKIND) :: snw_aging_table = 'snicar'

    ! debugging
    integer :: itmp1, itmp2
    real(kind=RKIND) :: rtmp1, rtmp2
    character(len=strKIND) :: ctmp1, ctmp2

    call MPAS_pool_get_config(domain % configs, "config_thermodynamics_type", config_thermodynamics_type)
    call MPAS_pool_get_config(domain % configs, "config_heat_conductivity_type", config_heat_conductivity_type)
    call MPAS_pool_get_config(domain % configs, "config_ocean_heat_transfer_type", config_ocean_heat_transfer_type)
    call MPAS_pool_get_config(domain % configs, "config_calc_surface_temperature", config_calc_surface_temperature)
    call MPAS_pool_get_config(domain % configs, "config_update_ocean_fluxes", config_update_ocean_fluxes)
    call MPAS_pool_get_config(domain % configs, "config_frazil_coupling_type", config_frazil_coupling_type)
    call MPAS_pool_get_config(domain % configs, "config_congelation_freezing_method", config_congelation_freezing_method)
    call MPAS_pool_get_config(domain % configs, "config_min_friction_velocity", config_min_friction_velocity)
    call MPAS_pool_get_config(domain % configs, "config_ice_ocean_drag_coefficient", config_ice_ocean_drag_coefficient)
    call MPAS_pool_get_config(domain % configs, "config_snow_thermal_conductivity", config_snow_thermal_conductivity)
    call MPAS_pool_get_config(domain % configs, "config_rapid_mode_channel_radius", config_rapid_mode_channel_radius)
    call MPAS_pool_get_config(domain % configs, "config_rapid_model_critical_Ra", config_rapid_model_critical_Ra)
    call MPAS_pool_get_config(domain % configs, "config_rapid_mode_aspect_ratio", config_rapid_mode_aspect_ratio)
    call MPAS_pool_get_config(domain % configs, "config_slow_mode_drainage_strength", config_slow_mode_drainage_strength)
    call MPAS_pool_get_config(domain % configs, "config_slow_mode_critical_porosity", config_slow_mode_critical_porosity)
    call MPAS_pool_get_config(domain % configs, "config_macro_drainage_timescale", config_macro_drainage_timescale)
    call MPAS_pool_get_config(domain % configs, "config_congelation_ice_porosity", config_congelation_ice_porosity)
!    call MPAS_pool_get_config(domain % configs, "config_frazil_ice_porosity", config_frazil_ice_porosity)
!    call MPAS_pool_get_config(domain % configs, "config_frazil_salinity_reduction", config_frazil_salinity_reduction)
    call MPAS_pool_get_config(domain % configs, "config_shortwave_type", config_shortwave_type)
    call MPAS_pool_get_config(domain % configs, "config_use_snicar_ad", config_use_snicar_ad)
    call MPAS_pool_get_config(domain % configs, "config_albedo_type", config_albedo_type)
    call MPAS_pool_get_config(domain % configs, "config_visible_ice_albedo", config_visible_ice_albedo)
    call MPAS_pool_get_config(domain % configs, "config_infrared_ice_albedo", config_infrared_ice_albedo)
    call MPAS_pool_get_config(domain % configs, "config_visible_snow_albedo", config_visible_snow_albedo)
    call MPAS_pool_get_config(domain % configs, "config_infrared_snow_albedo", config_infrared_snow_albedo)
    call MPAS_pool_get_config(domain % configs, "config_variable_albedo_thickness_limit", config_variable_albedo_thickness_limit)
!    call MPAS_pool_get_config(domain % configs, "config_snow_surface_scattering_layer_depth", config_snow_surface_scattering_layer_depth)
    call MPAS_pool_get_config(domain % configs, "config_ice_shortwave_tuning_parameter", config_ice_shortwave_tuning_parameter)
    call MPAS_pool_get_config(domain % configs, "config_pond_shortwave_tuning_parameter", config_pond_shortwave_tuning_parameter)
    call MPAS_pool_get_config(domain % configs, "config_snow_shortwave_tuning_parameter", config_snow_shortwave_tuning_parameter)
    call MPAS_pool_get_config(domain % configs, "config_use_shortwave_redistribution", config_use_shortwave_redistribution)
    call MPAS_pool_get_config(domain % configs, "config_shortwave_redistribution_fraction", config_shortwave_redistribution_fraction)
    call MPAS_pool_get_config(domain % configs, "config_shortwave_redistribution_threshold", config_shortwave_redistribution_threshold)
    call MPAS_pool_get_config(domain % configs, "config_temp_change_snow_grain_radius_change", &
                                                 config_temp_change_snow_grain_radius_change)
    call MPAS_pool_get_config(domain % configs, "config_max_melting_snow_grain_radius", config_max_melting_snow_grain_radius)
    call MPAS_pool_get_config(domain % configs, "config_algae_absorption_coefficient", config_algae_absorption_coefficient)
    call MPAS_pool_get_config(domain % configs, "config_ice_strength_formulation", config_ice_strength_formulation)
    call MPAS_pool_get_config(domain % configs, "config_ridging_participation_function", config_ridging_participation_function)
    call MPAS_pool_get_config(domain % configs, "config_ridging_redistribution_function", config_ridging_redistribution_function)
    call MPAS_pool_get_config(domain % configs, "config_ridging_efolding_scale", config_ridging_efolding_scale)
    call MPAS_pool_get_config(domain % configs, "config_ratio_ridging_work_to_PE", config_ratio_ridging_work_to_PE)
    call MPAS_pool_get_config(domain % configs, "config_atmos_boundary_method", config_atmos_boundary_method)
    call MPAS_pool_get_config(domain % configs, "config_calc_surface_stresses", config_calc_surface_stresses)
    call MPAS_pool_get_config(domain % configs, "config_use_form_drag", config_use_form_drag)
    call MPAS_pool_get_config(domain % configs, "config_use_high_frequency_coupling", config_use_high_frequency_coupling)
    call MPAS_pool_get_config(domain % configs, "config_boundary_layer_iteration_number", config_boundary_layer_iteration_number)
    call MPAS_pool_get_config(domain % configs, "config_use_ocean_mixed_layer", config_use_ocean_mixed_layer)
    call MPAS_pool_get_config(domain % configs, "config_sea_freezing_temperature_type", config_sea_freezing_temperature_type)
    call MPAS_pool_get_config(domain % configs, "config_itd_conversion_type", config_itd_conversion_type)
    call MPAS_pool_get_config(domain % configs, "config_category_bounds_type", config_category_bounds_type)
    call MPAS_pool_get_config(domain % configs, "config_snow_to_ice_transition_depth", config_snow_to_ice_transition_depth)
    call MPAS_pool_get_config(domain % configs, "config_pond_refreezing_type", config_pond_refreezing_type)
    call MPAS_pool_get_config(domain % configs, "config_salt_flux_coupling_type", config_salt_flux_coupling_type)
    call MPAS_pool_get_config(domain % configs, "config_pond_flushing_factor", config_pond_flushing_factor)
    call MPAS_pool_get_config(domain % configs, "config_min_meltwater_retained_fraction", config_min_meltwater_retained_fraction)
    call MPAS_pool_get_config(domain % configs, "config_max_meltwater_retained_fraction", config_max_meltwater_retained_fraction)
    call MPAS_pool_get_config(domain % configs, "config_pond_depth_to_fraction_ratio", config_pond_depth_to_fraction_ratio)
    call MPAS_pool_get_config(domain % configs, "config_snow_on_pond_ice_tapering_parameter", &
                                                 config_snow_on_pond_ice_tapering_parameter)
    call MPAS_pool_get_config(domain % configs, "config_critical_pond_ice_thickness", config_critical_pond_ice_thickness)
    call MPAS_pool_get_config(domain % configs, "config_skeletal_bgc_flux_type", config_skeletal_bgc_flux_type)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_tracers", config_use_vertical_tracers)
    call MPAS_pool_get_config(domain % configs, "config_scale_initial_vertical_bgc", config_scale_initial_vertical_bgc)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_shortwave_bioabsorption", config_use_shortwave_bioabsorption)
    call MPAS_pool_get_config(domain % configs, "config_use_modal_aerosols", config_use_modal_aerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_macromolecules", config_use_macromolecules)
    call MPAS_pool_get_config(domain % configs, "config_use_iron_solubility_file", config_use_iron_solubility_file)
    call MPAS_pool_get_config(domain % configs, "config_do_restart_bgc", config_do_restart_bgc)
    call MPAS_pool_get_config(domain % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_biogrid_bottom_molecular_sublayer", &
                                                 config_biogrid_bottom_molecular_sublayer)
    call MPAS_pool_get_config(domain % configs, "config_bio_gravity_drainage_length_scale", &
                                                 config_bio_gravity_drainage_length_scale)
    call MPAS_pool_get_config(domain % configs, "config_biogrid_top_molecular_sublayer", config_biogrid_top_molecular_sublayer)
    call MPAS_pool_get_config(domain % configs, "config_snow_porosity_at_ice_surface", config_snow_porosity_at_ice_surface)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Si_to_N_diatoms", config_ratio_Si_to_N_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Si_to_N_small_plankton", config_ratio_Si_to_N_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Si_to_N_phaeocystis", config_ratio_Si_to_N_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_ratio_S_to_N_diatoms", config_ratio_S_to_N_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_ratio_S_to_N_small_plankton", config_ratio_S_to_N_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_ratio_S_to_N_phaeocystis", config_ratio_S_to_N_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_C_diatoms", config_ratio_Fe_to_C_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_C_small_plankton", config_ratio_Fe_to_C_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_C_phaeocystis", config_ratio_Fe_to_C_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_N_diatoms", config_ratio_Fe_to_N_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_N_small_plankton", config_ratio_Fe_to_N_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_N_phaeocystis", config_ratio_Fe_to_N_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_DON", config_ratio_Fe_to_DON)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_DOC_saccharids", config_ratio_Fe_to_DOC_saccharids)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_DOC_lipids", config_ratio_Fe_to_DOC_lipids)
    call MPAS_pool_get_config(domain % configs, "config_respiration_fraction_of_growth", config_respiration_fraction_of_growth)
    call MPAS_pool_get_config(domain % configs, "config_algal_maximum_velocity", config_algal_maximum_velocity)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_dust", config_ratio_Fe_to_dust)
    call MPAS_pool_get_config(domain % configs, "config_solubility_of_Fe_in_dust", config_solubility_of_Fe_in_dust)
    call MPAS_pool_get_config(domain % configs, "config_chla_absorptivity_of_diatoms", config_chla_absorptivity_of_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_chla_absorptivity_of_small_plankton", &
                                                 config_chla_absorptivity_of_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_chla_absorptivity_of_phaeocystis", config_chla_absorptivity_of_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_light_attenuation_diatoms", config_light_attenuation_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_light_attenuation_small_plankton", config_light_attenuation_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_light_attenuation_phaeocystis", config_light_attenuation_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_light_inhibition_diatoms", config_light_inhibition_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_light_inhibition_small_plankton", config_light_inhibition_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_light_inhibition_phaeocystis", config_light_inhibition_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_maximum_growth_rate_diatoms", config_maximum_growth_rate_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_maximum_growth_rate_small_plankton", &
                                                 config_maximum_growth_rate_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_maximum_growth_rate_phaeocystis", config_maximum_growth_rate_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_temperature_growth_diatoms", config_temperature_growth_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_temperature_growth_small_plankton", &
                                                 config_temperature_growth_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_temperature_growth_phaeocystis", config_temperature_growth_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_grazed_fraction_diatoms", config_grazed_fraction_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_grazed_fraction_small_plankton", config_grazed_fraction_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_grazed_fraction_phaeocystis", config_grazed_fraction_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_mortality_diatoms", config_mortality_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_mortality_small_plankton", config_mortality_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_mortality_phaeocystis", config_mortality_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_temperature_mortality_diatoms", config_temperature_mortality_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_temperature_mortality_small_plankton", &
                                                 config_temperature_mortality_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_temperature_mortality_phaeocystis", &
                                                 config_temperature_mortality_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_exudation_diatoms", config_exudation_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_exudation_small_plankton", config_exudation_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_exudation_phaeocystis", config_exudation_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_nitrate_saturation_diatoms", config_nitrate_saturation_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_nitrate_saturation_small_plankton", &
                                                 config_nitrate_saturation_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_nitrate_saturation_phaeocystis", config_nitrate_saturation_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_ammonium_saturation_diatoms", config_ammonium_saturation_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_ammonium_saturation_small_plankton", &
                                                 config_ammonium_saturation_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_ammonium_saturation_phaeocystis", &
                                                 config_ammonium_saturation_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_silicate_saturation_diatoms", config_silicate_saturation_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_silicate_saturation_small_plankton", &
                                                 config_silicate_saturation_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_silicate_saturation_phaeocystis", config_silicate_saturation_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_iron_saturation_diatoms", config_iron_saturation_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_iron_saturation_small_plankton", config_iron_saturation_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_iron_saturation_phaeocystis", config_iron_saturation_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_fraction_spilled_to_DON", config_fraction_spilled_to_DON)
    call MPAS_pool_get_config(domain % configs, "config_degredation_of_DON", config_degredation_of_DON)
    call MPAS_pool_get_config(domain % configs, "config_fraction_DON_ammonium", config_fraction_DON_ammonium)
    call MPAS_pool_get_config(domain % configs, "config_fraction_loss_to_saccharids", config_fraction_loss_to_saccharids)
    call MPAS_pool_get_config(domain % configs, "config_fraction_loss_to_lipids",  config_fraction_loss_to_lipids)
    call MPAS_pool_get_config(domain % configs, "config_fraction_exudation_to_saccharids", config_fraction_exudation_to_saccharids)
    call MPAS_pool_get_config(domain % configs, "config_fraction_exudation_to_lipids", config_fraction_exudation_to_lipids)
    call MPAS_pool_get_config(domain % configs, "config_remineralization_saccharids", config_remineralization_saccharids)
    call MPAS_pool_get_config(domain % configs, "config_remineralization_lipids", config_remineralization_lipids)
    call MPAS_pool_get_config(domain % configs, "config_maximum_brine_temperature", config_maximum_brine_temperature)
    call MPAS_pool_get_config(domain % configs, "config_salinity_dependence_of_growth", config_salinity_dependence_of_growth)
    call MPAS_pool_get_config(domain % configs, "config_minimum_optical_depth", config_minimum_optical_depth)
    call MPAS_pool_get_config(domain % configs, "config_slopped_grazing_fraction", config_slopped_grazing_fraction)
    call MPAS_pool_get_config(domain % configs, "config_excreted_fraction", config_excreted_fraction)
    call MPAS_pool_get_config(domain % configs, "config_fraction_mortality_to_ammonium", config_fraction_mortality_to_ammonium)
    call MPAS_pool_get_config(domain % configs, "config_fraction_iron_remineralized", config_fraction_iron_remineralized)
    call MPAS_pool_get_config(domain % configs, "config_nitrification_rate", config_nitrification_rate)
    call MPAS_pool_get_config(domain % configs, "config_desorption_loss_particulate_iron", config_desorption_loss_particulate_iron)
    call MPAS_pool_get_config(domain % configs, "config_maximum_loss_fraction", config_maximum_loss_fraction)
    call MPAS_pool_get_config(domain % configs, "config_maximum_ratio_iron_to_saccharids", config_maximum_ratio_iron_to_saccharids)
    call MPAS_pool_get_config(domain % configs, "config_respiration_loss_to_DMSPd", config_respiration_loss_to_DMSPd)
    call MPAS_pool_get_config(domain % configs, "config_DMSP_to_DMS_conversion_fraction", config_DMSP_to_DMS_conversion_fraction)
    call MPAS_pool_get_config(domain % configs, "config_DMSP_to_DMS_conversion_time", config_DMSP_to_DMS_conversion_time)
    call MPAS_pool_get_config(domain % configs, "config_DMS_oxidation_time", config_DMS_oxidation_time)
    call MPAS_pool_get_config(domain % configs, "config_snow_redistribution_scheme", config_snow_redistribution_scheme)
    call MPAS_pool_get_config(domain % configs, "config_fallen_snow_radius", config_fallen_snow_radius)
    call MPAS_pool_get_config(domain % configs, "config_use_snow_liquid_ponds", config_use_snow_liquid_ponds)
    call MPAS_pool_get_config(domain % configs, "config_new_snow_density", config_new_snow_density)
    call MPAS_pool_get_config(domain % configs, "config_max_snow_density", config_max_snow_density)
    call MPAS_pool_get_config(domain % configs, "config_minimum_wind_compaction", config_minimum_wind_compaction)
    call MPAS_pool_get_config(domain % configs, "config_wind_compaction_factor", config_wind_compaction_factor)
    call MPAS_pool_get_config(domain % configs, "config_snow_redistribution_factor", config_snow_redistribution_factor)
    call MPAS_pool_get_config(domain % configs, "config_max_dry_snow_radius", config_max_dry_snow_radius)
    call MPAS_pool_get_config(domain % configs, "config_use_snow_grain_radius", config_use_snow_grain_radius)
    call MPAS_pool_get_config(domain % configs, "config_floediam", config_floediam)
    call MPAS_pool_get_config(domain % configs, "config_floeshape", config_floeshape)

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nGrainAgingTemperature", nGrainAgingTemperature)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nGrainAgingTempGradient", nGrainAgingTempGradient)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nGrainAgingSnowDensity", nGrainAgingSnowDensity)

    call MPAS_pool_get_subpool(domain % blocklist % structs, "snow", snow)
    call MPAS_pool_get_array(snow, "snowEmpiricalGrowthParameterTau", snowEmpiricalGrowthParameterTau)
    call MPAS_pool_get_array(snow, "snowEmpiricalGrowthParameterKappa", snowEmpiricalGrowthParameterKappa)
    call MPAS_pool_get_array(snow, "snowPropertyRate", snowPropertyRate)

    config_thermodynamics_type_int = config_cice_int("config_thermodynamics_type", config_thermodynamics_type)
    config_ice_strength_formulation_int = config_cice_int("config_ice_strength_formulation", config_ice_strength_formulation)
    config_ridging_participation_function_int = config_cice_int("config_ridging_participation_function", config_ridging_participation_function)
    config_ridging_redistribution_function_int = config_cice_int("config_ridging_redistribution_function", config_ridging_redistribution_function)
    config_itd_conversion_type_int = config_cice_int("config_itd_conversion_type", config_itd_conversion_type)
    config_category_bounds_type_int = config_cice_int("config_category_bounds_type", config_category_bounds_type)

    if (config_use_snicar_ad .and. (trim(config_shortwave_type(1:4)) == "dEdd")) then
       tmp_config_shortwave_type = 'dEdd_snicar_ad'
       tmp_config_snw_ssp_table  = 'snicar'
    else
       tmp_config_shortwave_type = config_shortwave_type
       tmp_config_snw_ssp_table  = 'test'
    endif

!echmod debugging
    call mpas_log_write(' ')
    call mpas_log_write(" ----- debugging parameters -----")
    call mpas_log_write(" ----- compare values prior to icepack init -----")

    rtmp1 = seaicePuny
    call icepack_query_parameters(puny_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('puny differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaicePi  ! defined above using pi from ice_constants_colpkg.F90
    call icepack_query_parameters(pi_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('pi differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceSecondsPerDay
    call icepack_query_parameters(secday_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('secday differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceDensitySnow
    call icepack_query_parameters(rhos_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('rhos differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceDensityIce
    call icepack_query_parameters(rhoi_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('rhoi differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceDensitySeaWater
    call icepack_query_parameters(rhow_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('rhow differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceAirSpecificHeat
    call icepack_query_parameters(cp_air_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('cp_air differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceIceSnowEmissivity
    call icepack_query_parameters(emissivity_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('emissivity differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceFreshIceSpecificHeat
    call icepack_query_parameters(cp_ice_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('cp_ice differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceSeaWaterSpecificHeat
    call icepack_query_parameters(cp_ocn_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('cp_ocn differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceIceSurfaceRoughness
    call icepack_query_parameters(iceruf_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('iceruf differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceOceanAlbedo
    call icepack_query_parameters(albocn_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('albocn differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceGravity
    call icepack_query_parameters(gravit_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('gravit differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceDensityFreshwater
    call icepack_query_parameters(rhofresh_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('rhofresh differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceZvir
    call icepack_query_parameters(zvir_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('zvir differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceVonKarmanConstant
    call icepack_query_parameters(vonkar_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('vonkar differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceWaterVaporSpecificHeat
    call icepack_query_parameters(cp_wv_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('cp_wv differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceStefanBoltzmann
    call icepack_query_parameters(stefan_boltzmann_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('stefan_boltzmann differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceReferenceSalinity
    call icepack_query_parameters(ice_ref_salinity_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('ice_ref_salinity differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceMaximumSalinity
    call icepack_query_parameters(saltmax_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('saltmax differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceMeltingTemperatureDepression
    call icepack_query_parameters(depressT_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('depressT differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceFreshWaterFreezingPoint
    call icepack_query_parameters(Tffresh_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('Tffresh differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceLatentHeatSublimation
    call icepack_query_parameters(Lsub_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('Lsub differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceLatentHeatVaporization
    call icepack_query_parameters(Lvap_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('Lvap differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceIceSurfaceMeltingTemperature
    call icepack_query_parameters(Timelt_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('Timelt differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceSnowSurfaceMeltingTemperature
    call icepack_query_parameters(Tsmelt_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('Tsmelt differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceiceStrengthConstantHiblerP
    call icepack_query_parameters(Pstar_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('Pstar differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceiceStrengthConstantHiblerC
    call icepack_query_parameters(Cstar_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('Cstar differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceStabilityReferenceHeight
    call icepack_query_parameters(zref_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('zref differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceSnowPatchiness
    call icepack_query_parameters(snowpatch_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('snowpatch differs $r $r',realArgs=(/rtmp1,rtmp2/))

    ctmp1 = config_frazil_coupling_type
    call icepack_query_parameters(cpl_frazil_out=ctmp2)
    if (trim(ctmp1) /= trim(ctmp2)) call mpas_log_write('cpl_frazil differs')

    ctmp1 = config_congelation_freezing_method
    call icepack_query_parameters(congel_freeze_out=ctmp2)
    if (trim(ctmp1) /= trim(ctmp2)) call mpas_log_write('congel_freeze differs')

!    rtmp1 = skeletalLayerThickness
!    call icepack_query_parameters(sk_l_out=rtmp2)
!    if (rtmp1 /= rtmp2) call mpas_log_write('sk_l differs $r $r',realArgs=(/rtmp1,rtmp2/))

    ! commented-out arguments take the default values from icepack_parameters.F90
    call icepack_init_parameters(&
         !argcheck_in             = , &               ! may not be needed in the driver
         puny_in                 = seaicePuny, &
         bignum_in               = seaiceBigNumber, &
         pi_in                   = seaicePi, &
         secday_in               = seaiceSecondsPerDay, &
         rhos_in                 = seaiceDensitySnow, &
         rhoi_in                 = seaiceDensityIce, &
         rhow_in                 = seaiceDensitySeaWater, &
         cp_air_in               = seaiceAirSpecificHeat, & !
         emissivity_in           = seaiceIceSnowEmissivity, & !
         cp_ice_in               = seaiceFreshIceSpecificHeat, & !
         cp_ocn_in               = seaiceSeaWaterSpecificHeat, &
         hfrazilmin_in           = seaiceFrazilMinimumThickness, &
         floediam_in             = config_floediam, &
         depressT_in             = seaiceMeltingTemperatureDepression, &
         dragio_in               = config_ice_ocean_drag_coefficient, & ! calc_dragio not implemented
         !thickness_ocn_layer1_in = , &        ! not yet implemented in MPAS-SI (stealth feature)
         !iceruf_ocn_in           = , &        ! under-ice roughness, not yet implemented
         albocn_in               = seaiceOceanAlbedo, &
         gravit_in               = seaiceGravity, &
         viscosity_dyn_in        = seaiceBrineDynamicViscosity, &
         Tocnfrz_in              = seaiceFreezingTemperatureConstant, &
         rhofresh_in             = seaiceDensityFreshwater, &
         zvir_in                 = seaiceZvir, &
         vonkar_in               = seaiceVonKarmanConstant, &
         cp_wv_in                = seaiceWaterVaporSpecificHeat, & !
         stefan_boltzmann_in     = seaiceStefanBoltzmann, &
         ice_ref_salinity_in     = seaiceReferenceSalinity, & !
         Tffresh_in              = seaiceFreshWaterFreezingPoint, &
         Lsub_in                 = seaiceLatentHeatSublimation, & !
         Lvap_in                 = seaiceLatentHeatVaporization, & !
         Timelt_in               = seaiceIceSurfaceMeltingTemperature, & !
         Tsmelt_in               = seaiceSnowSurfaceMeltingTemperature, & !
         iceruf_in               = seaiceIceSurfaceRoughness, & !
         Cf_in                   = config_ratio_ridging_work_to_PE, &
         Pstar_in                = seaiceIceStrengthConstantHiblerP, &
         Cstar_in                = seaiceIceStrengthConstantHiblerC, &
         kappav_in               = seaiceExtinctionCoef, &
         kice_in                 = seaiceFreshIceConductivity, &
         ksno_in                 = config_snow_thermal_conductivity, &
         zref_in                 = seaiceStabilityReferenceHeight, &
         hs_min_in               = seaiceSnowMinimumThickness, &
         snowpatch_in            = seaiceSnowPatchiness, & ! ccsm3 radiation scheme
         !rhosi_in                = , &        ! brine, zbgc, zsalinity
         sk_l_in                 = skeletalLayerThickness, & !
         R_gC2molC_in            = gramsCarbonPerMolCarbon, &
         saltmax_in              = seaiceMaximumSalinity, &
         phi_init_in             = seaiceFrazilIcePorosity, &
         !min_salin_in            = , &        ! ktherm=1, brine, zsalinity
         !salt_loss_in            = , &        ! brine, zsalinity
         !min_bgc_in              = , &        ! shortwave
         dSin0_frazil_in         = seaiceFrazilSalinityReduction, &
         hs_ssl_in               = seaiceSnowSurfaceScatteringLayer, &
         awtvdr_in               = seaiceAlbedoWtVisibleDirect, &
         awtidr_in               = seaiceAlbedoWtNearIRDirect, &
         awtvdf_in               = seaiceAlbedoWtVisibleDiffuse, &
         awtidf_in               = seaiceAlbedoWtNearIRDiffuse, &
         qqqice_in               = seaiceQsatQiceConstant, &
         TTTice_in               = seaiceQsatTiceConstant, &
         qqqocn_in               = seaiceQsatQocnConstant, &
         TTTocn_in               = seaiceQsatTocnConstant, &
         ktherm_in               = config_thermodynamics_type_int, &
         conduct_in              = config_heat_conductivity_type, &
         fbot_xfer_type_in       = config_ocean_heat_transfer_type, &
         calc_Tsfc_in            = config_calc_surface_temperature, &
         !dts_b_in                = , &        ! brine, zsalinity
         cpl_frazil_in           = config_frazil_coupling_type, &
         congel_freeze_in        = config_congelation_freezing_method, &
         update_ocn_f_in         = config_update_ocean_fluxes, &
         ustar_min_in            = config_min_friction_velocity, &
         a_rapid_mode_in         = config_rapid_mode_channel_radius, &
         Rac_rapid_mode_in       = config_rapid_model_critical_Ra, &
         aspect_rapid_mode_in    = config_rapid_mode_aspect_ratio, &
         dSdt_slow_mode_in       = config_slow_mode_drainage_strength, &
         phi_c_slow_mode_in      = config_slow_mode_critical_porosity, &
         phi_i_mushy_in          = config_congelation_ice_porosity, &
         tscale_pnd_drain_in     = config_macro_drainage_timescale, &
!         shortwave_in            = config_shortwave_type, &
         shortwave_in            = tmp_config_shortwave_type, &
         albedo_type_in          = config_albedo_type, &
         albsnowi_in             = config_infrared_snow_albedo, &
         albicev_in              = config_visible_ice_albedo, &
         albicei_in              = config_infrared_ice_albedo, &
         albsnowv_in             = config_visible_snow_albedo, &
         ahmax_in                = config_variable_albedo_thickness_limit, &
         R_ice_in                = config_ice_shortwave_tuning_parameter, &
         R_pnd_in                = config_pond_shortwave_tuning_parameter, &
         R_snw_in                = config_snow_shortwave_tuning_parameter, &
         dT_mlt_in               = config_temp_change_snow_grain_radius_change, &
         rsnw_mlt_in             = config_max_melting_snow_grain_radius, &
         kalg_in                 = config_algae_absorption_coefficient, &
         kstrength_in            = config_ice_strength_formulation_int, &
         krdg_partic_in          = config_ridging_participation_function_int, &
         krdg_redist_in          = config_ridging_redistribution_function_int, &
         mu_rdg_in               = config_ridging_efolding_scale, &
         atmbndy_in              = config_atmos_boundary_method, &
         calc_strair_in          = config_calc_surface_stresses, &
         formdrag_in             = config_use_form_drag, &
         highfreq_in             = config_use_high_frequency_coupling, &
         natmiter_in             = config_boundary_layer_iteration_number, &
         !atmiter_conv_in         = , &        ! not yet implemented in MPAS-SI (stealth feature)
         !calc_dragio_in          = , &        ! not yet implemented in MPAS-SI (stealth feature)
         tfrz_option_in          = config_sea_freezing_temperature_type, &
         kitd_in                 = config_itd_conversion_type_int, &
         kcatbound_in            = config_category_bounds_type_int, &
         hs0_in                  = config_snow_to_ice_transition_depth, &
         frzpnd_in               = config_pond_refreezing_type, &
         saltflux_option_in      = config_salt_flux_coupling_type, &
         floeshape_in            = config_floeshape, &
         !wave_spec_in            = , &        ! not yet implemented in MPAS-SI (stealth feature)
         !wave_spec_type_in       = , &        ! not yet implemented in MPAS-SI (stealth feature)
         !nfreq_in                = , &        ! not yet implemented in MPAS-SI (stealth feature)
         dpscale_in              = config_pond_flushing_factor, &
         rfracmin_in             = config_min_meltwater_retained_fraction, &
         rfracmax_in             = config_max_meltwater_retained_fraction, &
         pndaspect_in            = config_pond_depth_to_fraction_ratio, &
         hs1_in                  = config_snow_on_pond_ice_tapering_parameter, &
         hp1_in                  = config_critical_pond_ice_thickness, &
         bgc_flux_type_in        = config_skeletal_bgc_flux_type, &
         z_tracers_in            = config_use_vertical_tracers, &
         scale_bgc_in            = config_scale_initial_vertical_bgc, &
         solve_zbgc_in           = config_use_vertical_biochemistry, &
         dEdd_algae_in           = config_use_shortwave_bioabsorption, &
         modal_aero_in           = config_use_modal_aerosols, &
         use_macromolecules_in   = config_use_macromolecules, &
         restartbgc_in           = config_do_restart_bgc, &
         skl_bgc_in              = config_use_skeletal_biochemistry, &
         grid_o_in               = config_biogrid_bottom_molecular_sublayer, &
         l_sk_in                 = config_bio_gravity_drainage_length_scale, &
         grid_o_t_in             = config_biogrid_top_molecular_sublayer, &
         !grid_oS_in              = config_zsalinity_molecular_sublayer, &
         !l_skS_in                = config_zsalinity_gravity_drainage_scale, &
         phi_snow_in             = config_snow_porosity_at_ice_surface, &
         fr_resp_in              = config_respiration_fraction_of_growth, &
         algal_vel_in            = config_algal_maximum_velocity, &
         R_dFe2dust_in           = config_ratio_Fe_to_dust, &
         dustFe_sol_in           = config_solubility_of_Fe_in_dust, &
         use_atm_dust_iron_in    = .not.config_use_iron_solubility_file, &
         T_max_in                = config_maximum_brine_temperature, &
         fsal_in                 = config_salinity_dependence_of_growth, &
         op_dep_min_in           = config_minimum_optical_depth, &
         fr_graze_s_in           = config_slopped_grazing_fraction, &
         fr_graze_e_in           = config_excreted_fraction, &
         fr_mort2min_in          = config_fraction_mortality_to_ammonium, &
         fr_dFe_in               = config_fraction_iron_remineralized, &
         k_nitrif_in             = config_nitrification_rate, &
         t_iron_conv_in          = config_desorption_loss_particulate_iron, &
         max_loss_in             = config_maximum_loss_fraction, &
         max_dfe_doc1_in         = config_maximum_ratio_iron_to_saccharids, &
         fr_resp_s_in            = config_respiration_loss_to_DMSPd, &
         y_sk_DMS_in             = config_DMSP_to_DMS_conversion_fraction, &
         t_sk_conv_in            = config_DMSP_to_DMS_conversion_time, &
         t_sk_ox_in              = config_DMS_oxidation_time, &
         !conserv_check_in        = , &
         snwgrain_in             = config_use_snow_grain_radius, &
         snwredist_in            = config_snow_redistribution_scheme, &
         use_smliq_pnd_in        = config_use_snow_liquid_ponds, &
         rsnw_fall_in            = config_fallen_snow_radius, &
         rsnw_tmax_in            = config_max_dry_snow_radius, &
         rhosnew_in              = config_new_snow_density, &
         rhosmax_in              = config_max_snow_density, &
         windmin_in              = config_minimum_wind_compaction, &
         drhosdwind_in           = config_wind_compaction_factor, &
         snwlvlfac_in            = config_snow_redistribution_factor, &
         isnw_T_in               = nGrainAgingTemperature, &
         isnw_Tgrd_in            = nGrainAgingTempGradient, &
         isnw_rhos_in            = nGrainAgingSnowDensity, &
         snowage_tau_in          = snowEmpiricalGrowthParameterTau(:,:,:), &
         snowage_kappa_in        = snowEmpiricalGrowthParameterKappa(:,:,:), &
         snowage_drdt0_in        = snowPropertyRate(:,:,:), &
         snw_aging_table_in      = snw_aging_table, &
         snw_ssp_table_in        = tmp_config_snw_ssp_table &
         )

    call seaice_icepack_write_warnings(icepack_warnings_aborted())

    call mpas_log_write(" ----- compare values after icepack init -----")

    ! Lfresh is derived in Icepack, not sent in from driver
    rtmp1 = seaiceLatentHeatMelting
    call icepack_query_parameters(Lfresh_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('Lfresh differs $r $r',realArgs=(/rtmp1,rtmp2/))

    ! recheck values
    rtmp1 = seaiceAirSpecificHeat
    call icepack_query_parameters(cp_air_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('cp_air differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceIceSnowEmissivity
    call icepack_query_parameters(emissivity_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('emissivity differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceFreshIceSpecificHeat
    call icepack_query_parameters(cp_ice_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('cp_ice differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceSeaWaterSpecificHeat
    call icepack_query_parameters(cp_ocn_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('cp_ocn differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceZvir
    call icepack_query_parameters(zvir_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('zvir differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceLatentHeatSublimation
    call icepack_query_parameters(Lsub_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('Lsub differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = seaiceSnowPatchiness
    call icepack_query_parameters(snowpatch_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('snowpatch differs $r $r',realArgs=(/rtmp1,rtmp2/))

    ctmp1 = config_frazil_coupling_type
    call icepack_query_parameters(cpl_frazil_out=ctmp2)
    if (trim(ctmp1) /= trim(ctmp2)) call mpas_log_write('cpl_frazil differs')

    ctmp1 = config_congelation_freezing_method
    call icepack_query_parameters(congel_freeze_out=ctmp2)
    if (trim(ctmp1) /= trim(ctmp2)) call mpas_log_write('congel_freeze differs')

    call mpas_log_write(" ----- end debugging parameters -----")
    call mpas_log_write(' ')
!echmod debugging end

! parameters that have not yet been fully implemented in the list above
    !call icepack_init_parameters(& 
    !     shortwave             = config_shortwave_type, &         ! not working correctly above
    !     oceanmixed_ice        = config_use_ocean_mixed_layer, &  ! not used in Icepack/columnphysics (driver only)
    !-----------------------------------------------------------------------
    ! Parameters for thermodynamics
    !-----------------------------------------------------------------------

    ! ktherm:
    ! type of thermodynamics
    ! 1 = Bitz and Lipscomb 1999
    ! 2 = mushy layer theory
    !ktherm = config_cice_int("config_thermodynamics_type", config_thermodynamics_type)

    ! conduct:
    ! 'MU71' or 'bubbly'
    !conduct = config_heat_conductivity_type

    ! fbot_xfer_type:
    ! transfer coefficient type for ice-ocean heat flux
    !fbot_xfer_type = config_ocean_heat_transfer_type

    ! calc_Tsfc:
    ! if true, calculate surface temperature
    ! if false, Tsfc is computed elsewhere and
    ! atmos-ice fluxes are provided to CICE
    !calc_Tsfc = config_calc_surface_temperature

    ! ustar_min:
    ! minimum friction velocity for ice-ocean heat flux
    !ustar_min = config_min_friction_velocity

    ! mushy thermodynamics:

    ! a_rapid_mode:
    ! channel radius for rapid drainage mode (m)
    !a_rapid_mode = config_rapid_mode_channel_radius

    ! Rac_rapid_mode:
    ! critical rayleigh number for rapid drainage mode
    !Rac_rapid_mode = config_rapid_model_critical_Ra

    ! aspect_rapid_mode:
    ! aspect ratio for rapid drainage mode (larger=wider)
    !aspect_rapid_mode = config_rapid_mode_aspect_ratio

    ! dSdt_slow_mode:
    ! slow mode drainage strength (m s-1 K-1)
    !dSdt_slow_mode = config_slow_mode_drainage_strength

    ! phi_c_slow_mode:
    ! liquid fraction porosity cutoff for slow mode
    !phi_c_slow_mode = config_slow_mode_critical_porosity

    ! tscale_pnd_drain:
    ! Timescale for macroscopic drainage
    !tscale_pnd_drain = config_macro_drainage_timescale

    ! phi_i_mushy:
    ! liquid fraction of congelation ice
    !phi_i_mushy = config_congelation_ice_porosity

    ! phi_init:
    ! initial liquid fraction of frazil ice
    !phi_init = seaiceFrazilIcePorosity

    ! dSin0_frazil:
    ! initial frazil salinity reduction
    !dSin0_frazil = seaiceFrazilSalinityReduction

    ! congel_freeze:
    ! method used for congelation freezing
    !congel_freeze = config_congelation_freezing_method

    !-----------------------------------------------------------------------
    ! Parameters for radiation
    !-----------------------------------------------------------------------

    ! shortwave:
    ! shortwave method, 'default' ('ccsm3') or 'dEdd'
    !shortwave = config_shortwave_type

    ! albedo_type:
    ! albedo parameterization, 'default' ('ccsm3') or 'constant'
    ! shortwave='dEdd' overrides this parameter
    !albedo_type = config_albedo_type

    ! baseline albedos for ccsm3 shortwave, set in namelist

    ! albicev:
    ! visible ice albedo for h > ahmax
    !albicev = config_visible_ice_albedo

    ! albicei:
    ! near-ir ice albedo for h > ahmax
    !albicei = config_infrared_ice_albedo

    ! albsnowv:
    ! cold snow albedo, visible
    !albsnowv = config_visible_snow_albedo

    ! albsnowi:
    ! cold snow albedo, near IR
    !albsnowi = config_infrared_snow_albedo

    ! ahmax:
    ! thickness above which ice albedo is constant (m)
    !ahmax = config_variable_albedo_thickness_limit

    ! dEdd tuning parameters, set in namelist

    ! hs_ssl:
    ! snow surface scattering layer depth
    !hs_ssl = seaiceSnowSurfaceScatteringLayer

    ! R_ice:
    ! sea ice tuning parameter; +1 > 1sig increase in albedo
    !R_ice = config_ice_shortwave_tuning_parameter

    ! R_pnd:
    ! ponded ice tuning parameter; +1 > 1sig increase in albedo
    !R_pnd = config_pond_shortwave_tuning_parameter

    ! R_snw:
    ! snow tuning parameter; +1 > ~.01 change in broadband albedo
    !R_snw = config_snow_shortwave_tuning_parameter

    ! sw_redist
    ! Redistribute shortwave from layers near the melting temperature to the surface
    !sw_redist = config_use_shortwave_redistribution

    ! sw_frac
    ! Fraction of shortwave moved from layers near the melting temperature to the surface
    !sw_frac = config_shortwave_redistribution_fraction

    ! sw_dtemp
    ! Temperature threshold for moving shortwave to the surface
    !sw_dtemp = config_shortwave_redistribution_threshold

    ! dT_mlt:
    ! change in temp for non-melt to melt snow grain radius change (C)
    !dT_mlt = config_temp_change_snow_grain_radius_change

    ! rsnw_mlt:
    ! maximum melting snow grain radius (10^-6 m)
    !rsnw_mlt = config_max_melting_snow_grain_radius

    ! kalg:
    ! algae absorption coefficient for 0.5 m thick layer
    !kalg = config_algae_absorption_coefficient

    !-----------------------------------------------------------------------
    ! Parameters for ridging and strength
    !-----------------------------------------------------------------------

    ! kstrength:
    ! 0 for simple Hibler (1979) formulation
    ! 1 for Rothrock (1975) pressure formulation
    !kstrength = config_cice_int("config_ice_strength_formulation", config_ice_strength_formulation)

    ! krdg_partic:
    ! 0 for Thorndike et al. (1975) formulation
    ! 1 for exponential participation function
    !krdg_partic = config_cice_int("config_ridging_participation_function", config_ridging_participation_function)

    ! krdg_redist:
    ! 0 for Hibler (1980) formulation
    ! 1 for exponential redistribution function
    !krdg_redist = config_cice_int("config_ridging_redistribution_function", config_ridging_redistribution_function)

    ! mu_rdg:
    ! gives e-folding scale of ridged ice (m^.5)
    ! (krdg_redist = 1)
    !mu_rdg = config_ridging_efolding_scale

    ! Cf
    ! ratio of ridging work to PE change in ridging (kstrength = 1)
    !Cf = config_ratio_ridging_work_to_PE

    !-----------------------------------------------------------------------
    ! Parameters for atmosphere
    !-----------------------------------------------------------------------

    ! atmbndy:
    ! atmo boundary method, 'similarity' or 'constant' or 'mixed'
    !atmbndy = config_atmos_boundary_method

    ! calc_strair:
    ! if true, calculate wind stress components
    !calc_strair = config_calc_surface_stresses

    ! formdrag:
    ! if true, calculate form drag
    !formdrag = config_use_form_drag

    ! highfreq:
    ! if true, use high frequency coupling
    !highfreq = config_use_high_frequency_coupling

    ! natmiter:
    ! number of iterations for boundary layer calculations
    !natmiter = config_boundary_layer_iteration_number

    !-----------------------------------------------------------------------
    ! Parameters for ocean
    !-----------------------------------------------------------------------

    ! oceanmixed_ice:
    ! if true, use ocean mixed layer
    ! oceanmixed_ice = config_use_ocean_mixed_layer

    ! fbot_xfer_type:
    ! transfer coefficient type for ice-ocean heat flux
    ! fbot_xfer_type = config_ocean_heat_transfer_type

    ! tfrz_option:
    ! form of ocean freezing temperature
    ! 'minus1p8' = -1.8 C
    ! 'linear_salt' = -depressT * sss
    ! 'mushy' conforms with ktherm=2
    ! tfrz_option = config_sea_freezing_temperature_type

    ! dragio:
    ! Neutral ice-ocean drag coefficient
    ! dragio = config_ice_ocean_drag_coefficient

    !-----------------------------------------------------------------------
    ! Parameters for the ice thickness distribution
    !-----------------------------------------------------------------------

    ! kitd:
    ! type of itd conversions
    !   0 = delta function
    !   1 = linear remap
    !kitd = config_cice_int("config_itd_conversion_type", config_itd_conversion_type)

    ! kcatbound:
    !  -1 = single category
    !   0 = old category boundary formula
    !   1 = new formula giving round numbers
    !   2 = WMO standard
    !   3 = asymptotic formula
    !kcatbound = config_cice_int("config_category_bounds_type", config_category_bounds_type)

    !-----------------------------------------------------------------------
    ! Parameters for floe sizes
    !-----------------------------------------------------------------------

    ! floediam:
    ! effective floe diameter for lateral melt (m)
    ! default value is 300m from Steele (1992; doi:10.1029/92JC01755)
    !floediam = config_floediam


    ! floeshape:
    ! bulk floe shape constant for lateral melt (unitless coefficient)
    ! default value is 0.66 from Rothrock and Thorndike (1984; doi:10.1029/JC089iC04p06477)
    !floeshape = config_floeshape


    !-----------------------------------------------------------------------
    ! Parameters for melt ponds
    !-----------------------------------------------------------------------

    ! hs0:
    ! snow depth for transition to bare sea ice (m)
    !hs0 = config_snow_to_ice_transition_depth

    ! level-ice ponds

    ! frzpnd:
    ! pond refreezing parameterization
    !frzpnd = config_pond_refreezing_type

    ! dpscale:
    ! alters e-folding time scale for flushing with BL99 thermodynamics
    !dpscale = config_pond_flushing_factor

    ! rfracmin:
    ! minimum retained fraction of meltwater
    !rfracmin = config_min_meltwater_retained_fraction

    ! rfracmax:
    ! maximum retained fraction of meltwater
    !rfracmax = config_max_meltwater_retained_fraction

    ! pndaspect:
    ! ratio of pond depth to pond fraction
    !pndaspect = config_pond_depth_to_fraction_ratio

    ! hs1:
    ! tapering parameter for snow on pond ice
    !hs1 = config_snow_on_pond_ice_tapering_parameter

    ! topo ponds

    ! hp1
    ! critical parameter for pond ice thickness
    !hp1 = config_critical_pond_ice_thickness

    !-----------------------------------------------------------------------
    ! Parameters for biogeochemistry
    !-----------------------------------------------------------------------

    ! bgc_flux_type:
    ! bgc_flux_type = config_skeletal_bgc_flux_type

    ! z_tracers:
    ! if .true., bgc or aerosol tracers are vertically resolved
    !z_tracers = config_use_vertical_tracers

    ! scale_bgc:
    ! if .true., initialize bgc tracers proportionally with salinity
    !scale_bgc = config_scale_initial_vertical_bgc

    ! solve_zbgc:
    ! if .true., solve vertical biochemistry portion of code
    !solve_zbgc = config_use_vertical_biochemistry

    ! dEdd_algae:
    ! if .true., algal absorption of Shortwave is computed in the
    !dEdd_algae = config_use_shortwave_bioabsorption

    ! skl_bgc:
    ! if true, solve skeletal biochemistry
    !skl_bgc = config_use_skeletal_biochemistry

! zsalinity has been deprecated
!    ! solve_zsal:
!    ! if true, update salinity profile from solve_S_dt
!    !solve_zsal = config_use_vertical_zsalinity

    ! modal_aero:
    ! if true, use modal aerosal optical properties
    ! only for use with tr_aero or tr_zaero
    !modal_aero = config_use_shortwave_bioabsorption

    ! grid_o:
    ! for bottom flux
    !grid_o = config_biogrid_bottom_molecular_sublayer

    ! l_sk:
    ! characteristic diffusive scale (brine) (m)
    !l_sk =config_bio_gravity_drainage_length_scale

    ! grid_o_t:
    ! top grid point length scale
    !grid_o_t = config_biogrid_top_molecular_sublayer

    ! phi_snow:
    ! porosity of snow
    !phi_snow = config_snow_porosity_at_ice_surface

    ! initbio_frac:
    ! fraction of ocean tracer concentration used to initialize tracer
    !initbio_frac = config_new_ice_fraction_biotracer

    ! frazil_scav:
    ! multiple of ocean tracer concentration due to frazil scavenging
    !frazil_scav = config_fraction_biotracer_in_frazil

    ! ratio_Si2N_diatoms:
    ! ratio of algal Silicate to Nitrate (mol/mol)
    ! ratio_Si2N_diatoms = config_ratio_Si_to_N_diatoms

    ! ratio_Si2N_sp:
    ! ratio of algal Silicate to Nitrogen (mol/mol)
    ! ratio_Si2N_sp = config_ratio_Si_to_N_small_plankton

    ! ratio_Si2N_phaeo:
    ! ratio of algal Silicate to Nitrogen (mol/mol)
    ! ratio_Si2N_phaeo = config_ratio_Si_to_N_phaeocystis

    ! ratio_S2N_diatoms:
    ! ratio of algal Sulphur to Nitrogen (mol/mol)
    ! ratio_S2N_diatoms = config_ratio_S_to_N_diatoms

    ! ratio_S2N_sp:
    ! ratio of algal Sulphur to Nitrogen (mol/mol)
    ! ratio_S2N_sp = config_ratio_S_to_N_small_plankton

    ! ratio_S2N_phaeo:
    ! ratio of algal Sulphur to Nitrogen (mol/mol)
    ! ratio_S2N_phaeo = config_ratio_S_to_N_phaeocystis

    ! ratio_Fe2C_diatoms:
    ! ratio of algal iron to carbon (umol/mol)
    ! ratio_Fe2C_diatoms = config_ratio_Fe_to_C_diatoms

    ! ratio_Fe2C_sp:
    ! ratio of algal iron to carbon (umol/mol)
    ! ratio_Fe2C_sp = config_ratio_Fe_to_C_small_plankton

    ! ratio_Fe2C_phaeo:
    ! ratio of algal iron to carbon (umol/mol)
    ! ratio_Fe2C_phaeo = config_ratio_Fe_to_C_phaeocystis

    ! ratio_Fe2N_diatoms:
    ! ratio of algal iron to nitrogen (umol/mol)
    ! ratio_Fe2N_diatoms = config_ratio_Fe_to_N_diatoms

    ! ratio_Fe2N_sp:
    ! ratio of algal iron to nitrogen (umol/mol)
    ! ratio_Fe2N_sp = config_ratio_Fe_to_N_small_plankton

    ! ratio_Fe2N_phaeo:
    ! ratio of algal iron to nitrogen (umol/mol)
    ! ratio_Fe2N_phaeo = config_ratio_Fe_to_N_phaeocystis

    ! ratio_Fe2DON:
    ! ratio of iron to nitrogen of DON (nmol/umol)
    ! ratio_Fe2DON = config_ratio_Fe_to_DON

    ! ratio_Fe2DOC_s:
    ! ratio of iron to carbon of DOC (nmol/umol) saccharids
    ! ratio_Fe2DOC_s = config_ratio_Fe_to_DOC_saccharids

    ! ratio_Fe2DOC_l:
    ! ratio of iron to carbon of DOC (nmol/umol) lipids
    ! ratio_Fe2DOC_l = config_ratio_Fe_to_DOC_lipids

    ! fr_resp:
    ! fraction of algal growth lost due to respiration
    ! fr_resp = config_respiration_fraction_of_growth

    ! tau_min:
    ! rapid mobile to stationary exchanges (s) = 1.5 hours
    ! tau_min = config_rapid_mobile_to_stationary_time

    ! tau_max:
    ! long time mobile to stationary exchanges (s) = 2 days
    ! tau_max = config_long_mobile_to_stationary_time

    ! algal_vel:
    ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
    ! algal_vel = config_algal_maximum_velocity

    ! R_dFe2dust:
    !  g/g (3.5% content) Tagliabue 2009
    ! R_dFe2dust = config_ratio_Fe_to_dust

    ! dustFe_sol;
    ! solubility fraction
    ! dustFe_sol = config_solubility_of_Fe_in_dust

    ! use_atm_dust_sol;
    ! compute dust iron contribution in icepack
    ! use_atm_dust_sol = .not.config_use_iron_solubility_file

    ! chlabs_diatoms:
    ! chl absorption (1/m/(mg/m^3))
    ! chlabs_diatoms = config_chla_absorptivity_of_diatoms

    ! chlabs_sp:
    ! chl absorption (1/m/(mg/m^3))
    ! chlabs_sp = config_chla_absorptivity_of_small_plankton

    ! chlabs_phaeo:
    ! chl absorption (1/m/(mg/m^3))
    ! chlabs_phaeo = config_chla_absorptivity_of_phaeocystis

    ! alpha2max_low_diatoms:
    ! light limitation diatoms (1/(W/m^2))
    ! alpha2max_low_diatoms = config_light_attenuation_diatoms

    ! alpha2max_low_sp:
    ! light limitation small plankton (1/(W/m^2))
    ! alpha2max_low_sp = config_light_attenuation_small_plankton

    ! alpha2max_low_phaeo:
    ! light limitation phaeocystis (1/(W/m^2))
    ! alpha2max_low_phaeo = config_light_attenuation_phaeocystis

    ! beta2max_diatoms:
    ! light inhibition diatoms(1/(W/m^2))
    ! beta2max_diatoms = config_light_inhibition_diatoms

    ! beta2max_sp:
    ! light inhibition small plankton(1/(W/m^2))
    ! beta2max_sp = config_light_inhibition_small_plankton

    ! beta2max_phaeo:
    ! light inhibition phaeocystis (1/(W/m^2))
    ! beta2max_phaeo = config_light_inhibition_phaeocystis

    ! mu_max_diatoms:
    ! maximum growth rate diatoms (1/day)
    ! mu_max_diatoms = config_maximum_growth_rate_diatoms

    ! mu_max_sp:
    ! maximum growth rate small plankton (1/day)
    ! mu_max_sp = config_maximum_growth_rate_small plankton

    ! mu_max_phaeo:
    ! maximum growth rate phaeocystis (1/day)
    ! mu_max_phaeo = config_maximum_growth_rate_phaeocystis

    ! grow_Tdep_sp:
    ! Temperature dependence of growth small plankton (1/C)
    ! grow_Tdep_sp = config_temperature_growth_small_plankton

    ! grow_Tdep_phaeo:
    ! Temperature dependence of growth phaeocystis (1/C)
    ! grow_Tdep_phaeo = config_temperature_growth_phaeocystis

    ! fr_graze_diatoms:
    ! Fraction grazed diatoms
    ! fr_graze_diatoms = config_grazed_fraction_diatoms

    ! fr_graze_sp:
    ! Fraction grazed small_plankton
    ! fr_graze_sp = config_grazed_fraction_small_plankton

    ! fr_graze_phaeo:
    ! Fraction grazed phaeocystis
    ! fr_graze_phaeo = config_grazed_fraction_phaeocystis

    ! mort_pre_diatoms:
    ! Mortality diatoms (1/day)
    ! mort_pre_diatoms = config_mortality_diatoms

    ! mort_pre_sp:
    ! Mortality small_plankton (1/day)
    ! mort_pre_sp = config_mortality_small_plankton

    ! mort_pre_phaeo:
    ! Mortality phaeocystis (1/day)
    ! mort_pre_phaeo = config_mortality_phaeocystis

    ! mort_Tdep_diatoms:
    ! T dependence of mortality diatoms (1/C)
    ! mort_Tdep_diatoms = config_temperature_mortality_diatoms

    ! mort_Tdep_sp:
    ! T dependence of mortality small plankton (1/C)
    ! mort_Tdep_sp = config_temperature_mortality_small_plankton

    ! mort_Tdep_phaeo:
    ! T dependence of mortality phaeocystis (1/C)
    ! mort_Tdep_phaeo = config_temperature_mortality_phaeocystis

    ! k_exude_diatoms:
    ! algal exudation diatoms (1/d)
    ! k_exude_diatoms = config_exudation_diatoms

    ! k_exude_sp:
    ! algal exudation small_plankton (1/d)
    ! k_exude_sp = config_exudation_small_plankton

    ! k_exude_phaeo:
    ! algal exudation phaeocystis (1/d)
    ! k_exude_phaeo = config_exudation_phaeocystis

    ! K_Nit_diatoms:
    ! nitrate half saturation diatoms (mmol/m^3)
    ! K_Nit_diatoms = config_nitrate_saturation_diatoms

    ! K_Nit_sp:
    ! nitrate half saturation small_plankton (mmol/m^3)
    ! K_Nit_sp = config_nitrate_saturation_small_plankton

    ! K_Nit_phaeo:
    ! nitrate half saturation phaeocystis (mmol/m^3)
    ! K_Nit_phaeocystis = config_nitrate_saturation_phaeocystis

    ! K_Am_diatoms:
    ! ammonium half saturation diatoms (mmol/m^3)
    ! K_Am_diatoms = config_ammonium_saturation_diatoms

    ! K_Am_sp:
    ! ammonium half saturation small_plankton (mmol/m^3)
    ! K_Am_sp = config_ammonium_saturation_small_plankton

    ! K_Am_phaeo:
    ! ammonium half saturation phaeocystis (mmol/m^3)
    ! K_Am_phaeocystis = config_ammonium_saturation_phaeocystis

    ! K_Sil_diatoms:
    ! silicate half saturation diatoms (mmol/m^3)
    ! K_Sil_diatoms = config_silicate_saturation_diatoms

    ! K_Sil_sp:
    ! silicate half saturation small_plankton (mmol/m^3)
    ! K_Sil_sp = config_silicate_saturation_small_plankton

    ! K_Sil_phaeo:
    ! silicate half saturation phaeocystis (mmol/m^3)
    ! K_Sil_phaeocystis = config_silicate_saturation_phaeocystis

    ! K_Fe_diatoms:
    ! iron half saturation diatoms (nM)
    ! K_Fe_diatoms = config_iron_saturation_diatoms

    ! K_Fe_sp:
    ! iron half saturation small_plankton (nM)
    ! K_Fe_sp = config_iron_saturation_small_plankton

    ! K_Fe_phaeo:
    ! iron half saturation phaeocystis (nM)
    ! K_Fe_phaeocystis = config_iron_saturation_phaeocystis

    ! f_don_protein:
    ! fraction of spilled grazing to proteins    !
    ! f_don_protein = config_fraction_spilled_to_DON

    ! kn_bac_protein:
    ! Bacterial degredation of DON (1/d)    !     !
    ! kn_bac_protein = config_degredation_of_DON

    ! f_don_Am_protein:
    ! fraction of remineralized DON to ammonium    !
    ! f_don_Am_protein = config_fraction_DON_ammonium

    ! f_doc_s:
    ! fraction of mortality to DOC saccharids
    ! f_doc_s = config_fraction_loss_to_saccharids

    ! f_doc_l:
    ! fraction of mortality to DOC lipids
    ! f_doc_l = config_fraction_loss_to_lipids

    ! f_exude_s:
    ! fraction of exudation to DOC saccharids
    ! f_exude_s = config_fraction_exudation_to_saccharids

    ! f_exude_l:
    ! fraction of exudation to DOC lipids
    ! f_exude_l = config_fraction_exudation_to_lipids

    ! k_bac_s:
    ! Bacterial degredation of DOC (1/d) saccharids
    ! k_bac_s = config_remineralization_saccharids

    ! k_bac_l:
    ! Bacterial degredation of DOC (1/d) lipids
    ! k_bac_l = config_remineralization_lipids

    ! T_max:
    ! maximum temperature (C)
    ! T_max = config_maximum_brine_temperature

    ! fsal:
    ! Salinity limitation (ppt)
    ! fsal = config_salinity_dependence_of_growth

    ! op_dep_min:
    ! Light attenuates for optical depths exceeding min
    ! op_dep_min = config_minimum_optical_depth

    ! fr_graze_s:
    ! fraction of grazing spilled or slopped
    ! fr_graze_s = config_slopped_grazing_fraction

    ! fr_graze_e:
    ! fraction of assimilation excreted
    ! fr_graze_e = config_excreted_fraction

    ! fr_mort2min:
    ! fractionation of mortality to Am
    ! fr_mort2min = config_fraction_mortality_to_ammonium

    ! fr_dFe:
    ! remineralized nitrogen (in units of algal iron)
    ! fr_dFe = config_fraction_iron_remineralized

    ! k_nitrif:
    ! nitrification rate (1/day)
    ! k_nitrif = config_nitrification_rate

    ! t_iron_conv:
    ! desorption loss pFe to dFe (day)
    ! t_iron_conv = config_desorption_loss_particulate_iron

    ! max_loss:
    ! restrict uptake to % of remaining value
    ! max_loss = config_maximum_loss_fraction

    ! max_dfe_doc1:
    ! max ratio of dFe to saccharides in the ice  (nM Fe/muM C)
    ! max_dfe_doc1 = config_maximum_ratio_iron_to_saccharids

    ! fr_resp_s:
    ! DMSPd fraction of respiration loss as DMSPd
    ! fr_resp_s = config_respiration_loss_to_DMSPd

    ! y_sk_DMS:
    ! fraction conversion given high yield
    ! y_sk_DMS = config_DMSP_to_DMS_conversion_fraction

    ! t_sk_conv:
    ! Stefels conversion time (d)
    ! t_sk_conv = config_DMSP_to_DMS_conversion_time

    ! t_sk_ox:
    ! DMS oxidation time (d)
    ! t_sk_ox = config_DMS_oxidation_time

    ! algaltype_diatoms:
    ! mobility type diatoms
    ! algaltype_diatoms = config_mobility_type_diatoms

    ! algaltype_sp:
    ! mobility type small_plankton
    ! algaltype_sp = config_mobility_type_small_plankton

    ! algaltype_phaeo:
    ! mobility type phaeocystis
    ! algaltype_phaeo = config_mobility_type_phaeocystis

    ! nitratetype:
    ! mobility type nitrate
    ! nitratetype = config_mobility_type_nitrate

    ! ammoniumtype:
    ! mobility type ammonium
    ! ammoniumtype = config_mobility_type_ammonium

    ! silicatetype:
    ! mobility type silicate
    ! silicatetype = config_mobility_type_silicate

    ! dmspptype:
    ! mobility type DMSPp
    ! dmspptype = config_mobility_type_DMSPp

    ! dmspdtype:
    ! mobility type DMSPd
    ! dmspdtype = config_mobility_type_DMSPd

    ! humicstype:
    ! mobility type humics
    ! humicstype = config_mobility_type_humics

    ! doctype_s:
    ! mobility type sachharids
    ! doctype_s = config_mobility_type_saccharids

    ! doctype_l:
    ! mobility type lipids
    ! doctype_l = config_mobility_type_lipids

    ! dictype_1:
    ! mobility type dissolved inorganic carbon
    ! dictype_1 = config_mobility_type_inorganic_carbon

    ! dontype_protein:
    ! mobility type proteins
    ! dontype_protein = config_mobility_type_proteins

    ! fedtype_1:
    ! mobility type dissolved iron
    ! fedtype_1 =  config_mobility_type_dissolved_iron

    ! feptype_1:
    ! mobility type particulate iron
    ! feptype_1 =  config_mobility_type_particulate_iron

    ! zaerotype_bc1:
    ! mobility type for black carbon 1
    ! zaerotype_bc1 = config_mobility_type_black_carbon1

    ! zaerotype_bc2:
    ! mobility type for black carbon 2
    ! zaerotype_bc2 = config_mobility_type_black_carbon2

    ! zaerotype_dust1:
    ! mobility type for dust 1
    ! zaerotype_dust1 = config_mobility_type_dust1

    ! zaerotype_dust2:
    ! mobility type for dust 2
    ! zaerotype_dust2 = config_mobility_type_dust2

    ! zaerotype_dust3:
    ! mobility type for dust 3
    ! zaerotype_dust3 = config_mobility_type_dust3

    ! zaerotype_dust4:
    ! mobility type for dust 4
    ! zaerotype_dust4 = config_mobility_type_dust4

    ! ratio_C2N_diatoms:
    ! algal C to N ratio (mol/mol) diatoms
    ! ratio_C2N_diatoms = config_ratio_C_to_N_diatoms

    ! ratio_C2N_sp:
    ! algal C to N ratio (mol/mol) small_plankton
    ! ratio_C2N_sp = config_ratio_C_to_N_small_plankton

    ! ratio_C2N_phaeo:
    ! algal C to N ratio (mol/mol) phaeocystis
    ! ratio_C2N_phaeo = config_ratio_C_to_N_phaeocystis

    ! ratio_chl2N_diatoms:
    ! algal chla to N ratio (mol/mol) diatoms
    ! ratio_chl2N_diatoms = config_ratio_chla_to_N_diatoms

    ! ratio_chl2N_sp:
    ! algal chla to N ratio (mol/mol) small_plankton
    ! ratio_chl2N_sp = config_ratio_chla_to_N_small_plankton

    ! ratio_chl2N_phaeo:
    ! algal chla to N ratio (mol/mol) phaeocystis
    ! ratio_chl2N_phaeo = config_ratio_chla_to_N_phaeocystis

    ! F_abs_chl_diatoms:
    ! scales absorbed radiation for dEdd diatoms
    ! F_abs_chl_diatoms = config_scales_absorption_diatoms

    ! F_abs_chl_sp:
    ! scales absorbed radiation for dEdd small_plankton
    ! F_abs_chl_sp = config_scales_absorption_small_plankton

    ! F_abs_chl_phaeo:
    ! scales absorbed radiation for dEdd phaeocystis
    ! F_abs_chl_phaeo = config_scales_absorption_phaeocystis

    ! ratio_C2N_proteins:
    ! ratio of C to N in proteins (mol/mol)
    ! ratio_C2N_proteins = config_ratio_C_to_N_proteins

!    ! grid_oS:
!    ! for bottom flux (zsalinity)
!    !grid_oS = config_zsalinity_molecular_sublayer

!    ! l_skS:
!    ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)
!    !l_skS = config_zsalinity_gravity_drainage_scale

    !-----------------------------------------------------------------------
    ! Parameters for snow
    !-----------------------------------------------------------------------

    ! snwredist:
    ! snow redistribution type
    ! snwredist = config_snow_redistribution_scheme

    ! use_smliq_pnd:
    ! convert excess snow liquid to ponds
    ! use_smliq_pnd = config_use_snow_liquid_ponds

    ! rsnw_fall:
    ! fallen snow grain radius (um)
    ! rsnw_fall = config_fallen_snow_radius

    ! rsnw_tmax:
    ! maximum dry metamorphism snow grain radius (um)
    ! rsnw_tmax = config_max_dry_snow_radius

    ! rhosnew:
    ! new snow density (kg/m^3)
    ! rhosnew = config_new_snow_density

    ! rhosmax:
    ! maximum snow density (kg/m^3)
    ! rhosmax = config_max_snow_density

    ! windmin:
    ! minimum wind speed to compact snow (m/s)
    ! windmin = config_minimum_wind_compaction 

    ! drhosdwind:
    ! wind compaction factor (kg s/m^4)
    ! drhosdwind = config_wind_compaction_factor

    ! ksno:
    ! Thermal conductivity of snow (W/m/deg)
    ! ksno = config_snow_thermal_conductivity

  end subroutine init_icepack_package_configs

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  config_error
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 5th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine config_error(config_name, config_value, valid_options)

    character(len=*), intent(in) :: &
         config_name, &
         config_value, &
         valid_options

    call mpas_log_write("config_error: "//trim(config_name)//' has invalid value', messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(config_name)//': '//trim(config_value), messageType=MPAS_LOG_ERR)
    call mpas_log_write('valid options: '//trim(valid_options), messageType=MPAS_LOG_CRIT)

  end subroutine config_error

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  config_cice_int
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
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
       case ("BL99")
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

    ! kcatbound
    case ("config_category_bounds_type")

       select case (trim(configValue))
       case ("single category")
          configValueCice = -1
       case ("original")
          configValueCice = 0
       case ("new")
          configValueCice = 1
       case ("WMO")
          configValueCice = 2
       case ("asymptotic")
          configValueCice = 3
       end select

    ! kstrength
    case ("config_ice_strength_formulation")

       select case (trim(configValue))
       case ("Hibler79")
          configValueCice = 0
       case ("Rothrock75")
          configValueCice = 1
       end select

    ! krdg_partic
    case ("config_ridging_participation_function")

       select case (trim(configValue))
       case ("Thorndike75")
          configValueCice = 0
       case ("exponential")
          configValueCice = 1
       end select

    ! krdg_redist
    case ("config_ridging_redistribution_function")

       select case (trim(configValue))
       case ("Hibler80")
          configValueCice = 0
       case ("exponential")
          configValueCice = 1
       end select

    end select

  end function config_cice_int

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_icepack_non_activated_pointers
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_icepack_non_activated_pointers(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         drag, &
         tracers

    ! packages
    logical, pointer :: &
         pkgColumnTracerIceAgeActive, &
         pkgColumnTracerFirstYearIceActive, &
         pkgColumnTracerLevelIceActive, &
         pkgColumnTracerPondsActive, &
         pkgColumnTracerLidThicknessActive, &
         pkgColumnTracerAerosolsActive, &
         pkgColumnFormDragActive, &
         pkgColumnBiogeochemistryActive, &
         pkgTracerBrineActive, &
         pkgTracerMobileFractionActive, &
         pkgTracerSkeletalAlgaeActive, &
         pkgTracerSkeletalNitrateActive, &
         pkgTracerSkeletalCarbonActive, &
         pkgTracerSkeletalAmmoniumActive, &
         pkgTracerSkeletalSilicateActive, &
         pkgTracerSkeletalDMSActive, &
         pkgTracerSkeletalNonreactiveActive, &
         pkgTracerSkeletalHumicsActive, &
         pkgTracerSkeletalDONActive, &
         pkgTracerSkeletalIronActive, &
         pkgTracerVerticalAlgaeActive, &
         pkgTracerVerticalNitrateActive, &
         pkgTracerVerticalCarbonActive, &
         pkgTracerVerticalAmmoniumActive, &
         pkgTracerVerticalSilicateActive, &
         pkgTracerVerticalDMSActive, &
         pkgTracerVerticalNonreactiveActive, &
         pkgTracerVerticalHumicsActive, &
         pkgTracerVerticalDONActive, &
         pkgTracerVerticalIronActive, &
         pkgTracerZAerosolsActive, &
!         pkgTracerZSalinityActive, &           !echmod deprecate
         pkgColumnTracerEffectiveSnowDensityActive, &
         pkgColumnTracerSnowGrainRadiusActive


    ! mesh stand-ins
    type(field1DReal), pointer :: &
         latCell, lonCell ! nCells array

    type(field3DReal), pointer :: &
         iceAreaCategory

    ! drag variables
    type(field1DReal), pointer :: &
         oceanDragCoefficientSkin, &
         oceanDragCoefficientFloe, &
         oceanDragCoefficientKeel, &
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
         dragFloeSeparation

    block => domain % blocklist
    do while (associated(block))

       !-----------------------------------------------------------------------
       ! tracers
       !-----------------------------------------------------------------------

       call MPAS_pool_get_package(block % packages, "pkgColumnTracerIceAgeActive", pkgColumnTracerIceAgeActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerFirstYearIceActive", pkgColumnTracerFirstYearIceActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerLevelIceActive", pkgColumnTracerLevelIceActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerPondsActive", pkgColumnTracerPondsActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerLidThicknessActive", pkgColumnTracerLidThicknessActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerAerosolsActive", pkgColumnTracerAerosolsActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnBiogeochemistryActive", pkgColumnBiogeochemistryActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerBrineActive", pkgTracerBrineActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerMobileFractionActive", pkgTracerMobileFractionActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalAlgaeActive", pkgTracerSkeletalAlgaeActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalNitrateActive", pkgTracerSkeletalNitrateActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalCarbonActive", pkgTracerSkeletalCarbonActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalAmmoniumActive", pkgTracerSkeletalAmmoniumActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalSilicateActive", pkgTracerSkeletalSilicateActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalDMSActive", pkgTracerSkeletalDMSActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalNonreactiveActive", pkgTracerSkeletalNonreactiveActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalHumicsActive", pkgTracerSkeletalHumicsActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalDONActive", pkgTracerSkeletalDONActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalIronActive", pkgTracerSkeletalIronActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalAlgaeActive", pkgTracerVerticalAlgaeActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalNitrateActive", pkgTracerVerticalNitrateActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalCarbonActive", pkgTracerVerticalCarbonActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalAmmoniumActive", pkgTracerVerticalAmmoniumActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalSilicateActive", pkgTracerVerticalSilicateActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalDMSActive", pkgTracerVerticalDMSActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalNonreactiveActive", pkgTracerVerticalNonreactiveActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalHumicsActive", pkgTracerVerticalHumicsActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalDONActive", pkgTracerVerticalDONActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalIronActive", pkgTracerVerticalIronActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerZAerosolsActive", pkgTracerZAerosolsActive)
!       call MPAS_pool_get_package(block % packages, "pkgTracerZSalinityActive", pkgTracerZSalinityActive) !echmod - deprecate
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerEffectiveSnowDensityActive", pkgColumnTracerEffectiveSnowDensityActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerSnowGrainRadiusActive", pkgColumnTracerSnowGrainRadiusActive)


       ! ice age
       if (.not. pkgColumnTracerIceAgeActive) then
          call set_stand_in_tracer_array(block, "iceAge")
       endif

       ! first year ice
       if (.not. pkgColumnTracerFirstYearIceActive) then
          call set_stand_in_tracer_array(block, "firstYearIceArea")
       endif

       ! level ice
       if (.not. pkgColumnTracerLevelIceActive) then
          call set_stand_in_tracer_array(block, "levelIceArea")
          call set_stand_in_tracer_array(block, "levelIceVolume")
       endif

       ! ponds
       if (.not. pkgColumnTracerPondsActive) then
          call set_stand_in_tracer_array(block, "pondArea")
          call set_stand_in_tracer_array(block, "pondDepth")
       endif

       ! pond lids
       if (.not. pkgColumnTracerLidThicknessActive) then
          call set_stand_in_tracer_array(block, "pondLidThickness")
       endif

       ! aerosols
       if (.not. pkgColumnTracerAerosolsActive) then
          call set_stand_in_tracer_array(block, "snowScatteringAerosol")
          call set_stand_in_tracer_array(block, "snowBodyAerosol")
          call set_stand_in_tracer_array(block, "iceScatteringAerosol")
          call set_stand_in_tracer_array(block, "iceBodyAerosol")
       endif

       ! biogeochemistry
       if (.not. pkgTracerBrineActive) then
          call set_stand_in_tracer_array(block, "brineFraction")
       endif
       if (.not. pkgTracerMobileFractionActive) then
          call set_stand_in_tracer_array(block, "mobileFraction")
       endif
       if (.not. pkgTracerSkeletalAlgaeActive) then
          call set_stand_in_tracer_array(block, "skeletalAlgaeConc")
       endif
       if (.not. pkgTracerSkeletalNitrateActive) then
          call set_stand_in_tracer_array(block, "skeletalNitrateConc")
       endif
       if (.not. pkgTracerSkeletalSilicateActive) then
          call set_stand_in_tracer_array(block, "skeletalSilicateConc")
       endif
       if (.not. pkgTracerSkeletalAmmoniumActive) then
             call set_stand_in_tracer_array(block, "skeletalAmmoniumConc")
       endif
       if (.not. pkgTracerSkeletalDMSActive) then
          call set_stand_in_tracer_array(block, "skeletalDMSPpConc")
          call set_stand_in_tracer_array(block, "skeletalDMSPpConc")
          call set_stand_in_tracer_array(block, "skeletalDMSConc")
       endif
       if (.not. pkgTracerSkeletalCarbonActive) then
          call set_stand_in_tracer_array(block, "skeletalDOCConc")
          call set_stand_in_tracer_array(block, "skeletalDICConc")
       endif
       if (.not. pkgTracerSkeletalDONActive) then
          call set_stand_in_tracer_array(block, "skeletalDONConc")
       endif
       if (.not. pkgTracerSkeletalNonreactiveActive) then
          call set_stand_in_tracer_array(block, "skeletalNonreactiveConc")
       endif
       if (.not. pkgTracerSkeletalHumicsActive) then
          call set_stand_in_tracer_array(block, "skeletalHumicsConc")
       endif
       if (.not. pkgTracerSkeletalIronActive) then
          call set_stand_in_tracer_array(block, "skeletalParticulateIronConc")
          call set_stand_in_tracer_array(block, "skeletalDissolvedIronConc")
       endif
       if (.not. pkgTracerVerticalAlgaeActive) then
          call set_stand_in_tracer_array(block, "verticalAlgaeConc")
          call set_stand_in_tracer_array(block, "verticalAlgaeSnow")
          call set_stand_in_tracer_array(block, "verticalAlgaeIce")
       endif
       if (.not. pkgTracerVerticalNitrateActive) then
          call set_stand_in_tracer_array(block, "verticalNitrateConc")
          call set_stand_in_tracer_array(block, "verticalNitrateSnow")
          call set_stand_in_tracer_array(block, "verticalNitrateIce")
       endif
       if (.not. pkgTracerVerticalSilicateActive) then
          call set_stand_in_tracer_array(block, "verticalSilicateConc")
          call set_stand_in_tracer_array(block, "verticalSilicateSnow")
          call set_stand_in_tracer_array(block, "verticalSilicateIce")
       endif
       if (.not. pkgTracerVerticalAmmoniumActive) then
          call set_stand_in_tracer_array(block, "verticalAmmoniumConc")
          call set_stand_in_tracer_array(block, "verticalAmmoniumSnow")
          call set_stand_in_tracer_array(block, "verticalAmmoniumIce")
       endif
       if (.not. pkgTracerVerticalDMSActive) then
          call set_stand_in_tracer_array(block, "verticalDMSPpConc")
          call set_stand_in_tracer_array(block, "verticalDMSPdConc")
          call set_stand_in_tracer_array(block, "verticalDMSConc")
          call set_stand_in_tracer_array(block, "verticalDMSPpSnow")
          call set_stand_in_tracer_array(block, "verticalDMSPdSnow")
          call set_stand_in_tracer_array(block, "verticalDMSSnow")
          call set_stand_in_tracer_array(block, "verticalDMSPpIce")
          call set_stand_in_tracer_array(block, "verticalDMSPdIce")
          call set_stand_in_tracer_array(block, "verticalDMSIce")
       endif
       if (.not. pkgTracerVerticalCarbonActive) then
          call set_stand_in_tracer_array(block, "verticalDOCConc")
          call set_stand_in_tracer_array(block, "verticalDICConc")
          call set_stand_in_tracer_array(block, "verticalDOCSnow")
          call set_stand_in_tracer_array(block, "verticalDICSnow")
          call set_stand_in_tracer_array(block, "verticalDOCIce")
          call set_stand_in_tracer_array(block, "verticalDICIce")
       endif
       if (.not. pkgTracerVerticalDONActive) then
          call set_stand_in_tracer_array(block, "verticalDONConc")
          call set_stand_in_tracer_array(block, "verticalDONSnow")
          call set_stand_in_tracer_array(block, "verticalDONIce")
       endif
       if (.not. pkgTracerVerticalNonreactiveActive) then
          call set_stand_in_tracer_array(block, "verticalNonreactiveConc")
          call set_stand_in_tracer_array(block, "verticalNonreactiveSnow")
          call set_stand_in_tracer_array(block, "verticalNonreactiveIce")
       endif
       if (.not. pkgTracerVerticalHumicsActive) then
          call set_stand_in_tracer_array(block, "verticalHumicsConc")
          call set_stand_in_tracer_array(block, "verticalHumicsSnow")
          call set_stand_in_tracer_array(block, "verticalHumicsIce")
       endif
       if (.not. pkgTracerVerticalIronActive) then
          call set_stand_in_tracer_array(block, "verticalParticulateIronConc")
          call set_stand_in_tracer_array(block, "verticalDissolvedIronConc")
          call set_stand_in_tracer_array(block, "verticalParticulateIronSnow")
          call set_stand_in_tracer_array(block, "verticalDissolvedIronSnow")
          call set_stand_in_tracer_array(block, "verticalParticulateIronIce")
          call set_stand_in_tracer_array(block, "verticalDissolvedIronIce")
       endif
       if (.not. pkgTracerZAerosolsActive) then
          call set_stand_in_tracer_array(block, "verticalAerosolsConc")
          call set_stand_in_tracer_array(block, "verticalAerosolsSnow")
          call set_stand_in_tracer_array(block, "verticalAerosolsIce")
       endif
!       if (.not. pkgTracerZSalinityActive) then                          !echmod - deprecate
!          call set_stand_in_tracer_array(block, "verticalSalinity")
!       endif

       ! snow density tracer
       if (.not. pkgColumnTracerEffectiveSnowDensityActive) then
          call set_stand_in_tracer_array(block, "snowDensity")
       endif

       ! snow grain radius
       if (.not. pkgColumnTracerSnowGrainRadiusActive) then
          call set_stand_in_tracer_array(block, "snowLiquidMass")
          call set_stand_in_tracer_array(block, "snowIceMass")
          call set_stand_in_tracer_array(block, "snowGrainRadius")
       endif
       !-----------------------------------------------------------------------
       ! other column packages
       !-----------------------------------------------------------------------

       ! form drag
       call MPAS_pool_get_package(block % packages, "pkgColumnFormDragActive", pkgColumnFormDragActive)

       if (.not. pkgColumnFormDragActive) then

          ! get mesh stand-ins if have ones of right dimensions
          call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
          call MPAS_pool_get_field(mesh, "latCell", latCell) ! nCells real array

          call MPAS_pool_get_subpool(block % structs, "drag", drag)

          call MPAS_pool_get_field(drag, "oceanDragCoefficientSkin", oceanDragCoefficientSkin)
          call MPAS_pool_get_field(drag, "oceanDragCoefficientFloe", oceanDragCoefficientFloe)
          call MPAS_pool_get_field(drag, "oceanDragCoefficientKeel", oceanDragCoefficientKeel)
          call MPAS_pool_get_field(drag, "airDragCoefficientSkin", airDragCoefficientSkin)
          call MPAS_pool_get_field(drag, "airDragCoefficientFloe", airDragCoefficientFloe)
          call MPAS_pool_get_field(drag, "airDragCoefficientPond", airDragCoefficientPond)
          call MPAS_pool_get_field(drag, "airDragCoefficientRidge", airDragCoefficientRidge)
          call MPAS_pool_get_field(drag, "dragFreeboard", dragFreeboard)
          call MPAS_pool_get_field(drag, "dragIceSnowDraft", dragIceSnowDraft)
          call MPAS_pool_get_field(drag, "dragRidgeHeight", dragRidgeHeight)
          call MPAS_pool_get_field(drag, "dragRidgeSeparation", dragRidgeSeparation)
          call MPAS_pool_get_field(drag, "dragKeelDepth", dragKeelDepth)
          call MPAS_pool_get_field(drag, "dragKeelSeparation", dragKeelSeparation)
          call MPAS_pool_get_field(drag, "dragFloeLength", dragFloeLength)
          call MPAS_pool_get_field(drag, "dragFloeSeparation", dragFloeSeparation)

          oceanDragCoefficientSkin % array     => latCell % array
          oceanDragCoefficientFloe % array     => latCell % array
          oceanDragCoefficientKeel % array     => latCell % array
          airDragCoefficientSkin % array       => latCell % array
          airDragCoefficientFloe % array       => latCell % array
          airDragCoefficientPond % array       => latCell % array
          airDragCoefficientRidge % array      => latCell % array
          dragFreeboard % array                => latCell % array
          dragIceSnowDraft % array             => latCell % array
          dragRidgeHeight % array              => latCell % array
          dragRidgeSeparation % array          => latCell % array
          dragKeelDepth % array                => latCell % array
          dragKeelSeparation % array           => latCell % array
          dragFloeLength % array               => latCell % array
          dragFloeSeparation % array           => latCell % array

       endif

       block => block % next
    end do

  end subroutine init_icepack_non_activated_pointers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  finalize_icepack_non_activated_pointers
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine finalize_icepack_non_activated_pointers(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         drag, &
         tracers

    ! packages
    logical, pointer :: &
         pkgColumnTracerIceAgeActive, &
         pkgColumnTracerFirstYearIceActive, &
         pkgColumnTracerLevelIceActive, &
         pkgColumnTracerPondsActive, &
         pkgColumnTracerLidThicknessActive, &
         pkgColumnTracerAerosolsActive, &
         pkgColumnFormDragActive, &
         pkgColumnBiogeochemistryActive, &
         pkgTracerBrineActive, &
         pkgTracerMobileFractionActive, &
         pkgTracerSkeletalAlgaeActive, &
         pkgTracerSkeletalNitrateActive, &
         pkgTracerSkeletalCarbonActive, &
         pkgTracerSkeletalAmmoniumActive, &
         pkgTracerSkeletalSilicateActive, &
         pkgTracerSkeletalDMSActive, &
         pkgTracerSkeletalNonreactiveActive, &
         pkgTracerSkeletalHumicsActive, &
         pkgTracerSkeletalDONActive, &
         pkgTracerSkeletalIronActive, &
         pkgTracerVerticalAlgaeActive, &
         pkgTracerVerticalNitrateActive, &
         pkgTracerVerticalCarbonActive, &
         pkgTracerVerticalAmmoniumActive, &
         pkgTracerVerticalSilicateActive, &
         pkgTracerVerticalDMSActive, &
         pkgTracerVerticalNonreactiveActive, &
         pkgTracerVerticalHumicsActive, &
         pkgTracerVerticalDONActive, &
         pkgTracerVerticalIronActive, &
         pkgTracerZAerosolsActive, &
!         pkgTracerZSalinityActive, &                   !echmod deprecate
         pkgColumnTracerEffectiveSnowDensityActive, &
         pkgColumnTracerSnowGrainRadiusActive

    ! drag variables
    type(field1DReal), pointer :: &
         oceanDragCoefficientSkin, &
         oceanDragCoefficientFloe, &
         oceanDragCoefficientKeel, &
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
         dragFloeSeparation

    block => domain % blocklist
    do while (associated(block))

       !-----------------------------------------------------------------------
       ! tracers
       !-----------------------------------------------------------------------

       call MPAS_pool_get_package(block % packages, "pkgColumnTracerIceAgeActive", pkgColumnTracerIceAgeActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerFirstYearIceActive", pkgColumnTracerFirstYearIceActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerLevelIceActive", pkgColumnTracerLevelIceActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerPondsActive", pkgColumnTracerPondsActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerLidThicknessActive", pkgColumnTracerLidThicknessActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerAerosolsActive", pkgColumnTracerAerosolsActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnBiogeochemistryActive", pkgColumnBiogeochemistryActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerBrineActive", pkgTracerBrineActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerMobileFractionActive", pkgTracerMobileFractionActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalAlgaeActive", pkgTracerSkeletalAlgaeActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalNitrateActive", pkgTracerSkeletalNitrateActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalCarbonActive", pkgTracerSkeletalCarbonActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalAmmoniumActive", pkgTracerSkeletalAmmoniumActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalSilicateActive", pkgTracerSkeletalSilicateActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalDMSActive", pkgTracerSkeletalDMSActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalNonreactiveActive", pkgTracerSkeletalNonreactiveActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalHumicsActive", pkgTracerSkeletalHumicsActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalDONActive", pkgTracerSkeletalDONActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerSkeletalIronActive", pkgTracerSkeletalIronActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalAlgaeActive", pkgTracerVerticalAlgaeActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalNitrateActive", pkgTracerVerticalNitrateActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalCarbonActive", pkgTracerVerticalCarbonActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalAmmoniumActive", pkgTracerVerticalAmmoniumActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalSilicateActive", pkgTracerVerticalSilicateActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalDMSActive", pkgTracerVerticalDMSActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalNonreactiveActive", pkgTracerVerticalNonreactiveActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalHumicsActive", pkgTracerVerticalHumicsActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalDONActive", pkgTracerVerticalDONActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerVerticalIronActive", pkgTracerVerticalIronActive)
       call MPAS_pool_get_package(block % packages, "pkgTracerZAerosolsActive", pkgTracerZAerosolsActive)
!       call MPAS_pool_get_package(block % packages, "pkgTracerZSalinityActive", pkgTracerZSalinityActive)      !echmod deprecate
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerEffectiveSnowDensityActive", pkgColumnTracerEffectiveSnowDensityActive)
       call MPAS_pool_get_package(block % packages, "pkgColumnTracerSnowGrainRadiusActive", pkgColumnTracerSnowGrainRadiusActive)

       ! ice age
       if (.not. pkgColumnTracerIceAgeActive) then
          call finalize_stand_in_tracer_array(block, "iceAge")
       endif

       ! first year ice
       if (.not. pkgColumnTracerFirstYearIceActive) then
          call finalize_stand_in_tracer_array(block, "firstYearIceArea")
       endif

       ! level ice
       if (.not. pkgColumnTracerLevelIceActive) then
          call finalize_stand_in_tracer_array(block, "levelIceArea")
          call finalize_stand_in_tracer_array(block, "levelIceVolume")
       endif

       ! ponds
       if (.not. pkgColumnTracerPondsActive) then
          call finalize_stand_in_tracer_array(block, "pondArea")
          call finalize_stand_in_tracer_array(block, "pondDepth")
       endif

       ! pond lids
       if (.not. pkgColumnTracerLidThicknessActive) then
          call finalize_stand_in_tracer_array(block, "pondLidThickness")
       endif

       ! aerosols
       if (.not. pkgColumnTracerAerosolsActive) then
          call finalize_stand_in_tracer_array(block, "snowScatteringAerosol")
          call finalize_stand_in_tracer_array(block, "snowBodyAerosol")
          call finalize_stand_in_tracer_array(block, "iceScatteringAerosol")
          call finalize_stand_in_tracer_array(block, "iceBodyAerosol")
       endif

       ! biogeochemistry
       if (.not. pkgTracerBrineActive) then
          call finalize_stand_in_tracer_array(block, "brineFraction")
       endif
       if (.not. pkgTracerMobileFractionActive) then
          call finalize_stand_in_tracer_array(block, "mobileFraction")
       endif
       if (.not. pkgTracerSkeletalAlgaeActive) then
          call finalize_stand_in_tracer_array(block, "skeletalAlgaeConc")
       endif
       if (.not. pkgTracerSkeletalNitrateActive) then
          call finalize_stand_in_tracer_array(block, "skeletalNitrateConc")
       endif
       if (.not. pkgTracerSkeletalSilicateActive) then
          call finalize_stand_in_tracer_array(block, "skeletalSilicateConc")
       endif
       if (.not. pkgTracerSkeletalAmmoniumActive) then
             call finalize_stand_in_tracer_array(block, "skeletalAmmoniumConc")
       endif
       if (.not. pkgTracerSkeletalDMSActive) then
          call finalize_stand_in_tracer_array(block, "skeletalDMSPpConc")
          call finalize_stand_in_tracer_array(block, "skeletalDMSPpConc")
          call finalize_stand_in_tracer_array(block, "skeletalDMSConc")
       endif
       if (.not. pkgTracerSkeletalCarbonActive) then
          call finalize_stand_in_tracer_array(block, "skeletalDOCConc")
          call finalize_stand_in_tracer_array(block, "skeletalDICConc")
       endif
       if (.not. pkgTracerSkeletalDONActive) then
          call finalize_stand_in_tracer_array(block, "skeletalDONConc")
       endif
       if (.not. pkgTracerSkeletalNonreactiveActive) then
          call finalize_stand_in_tracer_array(block, "skeletalNonreactiveConc")
       endif
       if (.not. pkgTracerSkeletalHumicsActive) then
          call finalize_stand_in_tracer_array(block, "skeletalHumicsConc")
       endif
       if (.not. pkgTracerSkeletalIronActive) then
          call finalize_stand_in_tracer_array(block, "skeletalParticulateIronConc")
          call finalize_stand_in_tracer_array(block, "skeletalDissolvedIronConc")
       endif
       if (.not. pkgTracerVerticalAlgaeActive) then
          call finalize_stand_in_tracer_array(block, "verticalAlgaeConc")
          call finalize_stand_in_tracer_array(block, "verticalAlgaeSnow")
          call finalize_stand_in_tracer_array(block, "verticalAlgaeIce")
       endif
       if (.not. pkgTracerVerticalNitrateActive) then
          call finalize_stand_in_tracer_array(block, "verticalNitrateConc")
          call finalize_stand_in_tracer_array(block, "verticalNitrateSnow")
          call finalize_stand_in_tracer_array(block, "verticalNitrateIce")
       endif
       if (.not. pkgTracerVerticalSilicateActive) then
          call finalize_stand_in_tracer_array(block, "verticalSilicateConc")
          call finalize_stand_in_tracer_array(block, "verticalSilicateSnow")
          call finalize_stand_in_tracer_array(block, "verticalSilicateIce")
       endif
       if (.not. pkgTracerVerticalAmmoniumActive) then
          call finalize_stand_in_tracer_array(block, "verticalAmmoniumConc")
          call finalize_stand_in_tracer_array(block, "verticalAmmoniumSnow")
          call finalize_stand_in_tracer_array(block, "verticalAmmoniumIce")
       endif
       if (.not. pkgTracerVerticalDMSActive) then
          call finalize_stand_in_tracer_array(block, "verticalDMSPpConc")
          call finalize_stand_in_tracer_array(block, "verticalDMSPdConc")
          call finalize_stand_in_tracer_array(block, "verticalDMSConc")
          call finalize_stand_in_tracer_array(block, "verticalDMSPpSnow")
          call finalize_stand_in_tracer_array(block, "verticalDMSPdSnow")
          call finalize_stand_in_tracer_array(block, "verticalDMSSnow")
          call finalize_stand_in_tracer_array(block, "verticalDMSPpIce")
          call finalize_stand_in_tracer_array(block, "verticalDMSPdIce")
          call finalize_stand_in_tracer_array(block, "verticalDMSIce")
       endif
       if (.not. pkgTracerVerticalCarbonActive) then
          call finalize_stand_in_tracer_array(block, "verticalDOCConc")
          call finalize_stand_in_tracer_array(block, "verticalDICConc")
          call finalize_stand_in_tracer_array(block, "verticalDOCSnow")
          call finalize_stand_in_tracer_array(block, "verticalDICSnow")
          call finalize_stand_in_tracer_array(block, "verticalDOCIce")
          call finalize_stand_in_tracer_array(block, "verticalDICIce")
       endif
       if (.not. pkgTracerVerticalDONActive) then
          call finalize_stand_in_tracer_array(block, "verticalDONConc")
          call finalize_stand_in_tracer_array(block, "verticalDONSnow")
          call finalize_stand_in_tracer_array(block, "verticalDONIce")
       endif
       if (.not. pkgTracerVerticalNonreactiveActive) then
          call finalize_stand_in_tracer_array(block, "verticalNonreactiveConc")
          call finalize_stand_in_tracer_array(block, "verticalNonreactiveSnow")
          call finalize_stand_in_tracer_array(block, "verticalNonreactiveIce")
       endif
       if (.not. pkgTracerVerticalHumicsActive) then
          call finalize_stand_in_tracer_array(block, "verticalHumicsConc")
          call finalize_stand_in_tracer_array(block, "verticalHumicsSnow")
          call finalize_stand_in_tracer_array(block, "verticalHumicsIce")
       endif
       if (.not. pkgTracerVerticalIronActive) then
          call finalize_stand_in_tracer_array(block, "verticalParticulateIronConc")
          call finalize_stand_in_tracer_array(block, "verticalDissolvedIronConc")
          call finalize_stand_in_tracer_array(block, "verticalParticulateIronSnow")
          call finalize_stand_in_tracer_array(block, "verticalDissolvedIronSnow")
          call finalize_stand_in_tracer_array(block, "verticalParticulateIronIce")
          call finalize_stand_in_tracer_array(block, "verticalDissolvedIronIce")
       endif
       if (.not. pkgTracerZAerosolsActive) then
          call finalize_stand_in_tracer_array(block, "verticalAerosolsConc")
          call finalize_stand_in_tracer_array(block, "verticalAerosolsSnow")
          call finalize_stand_in_tracer_array(block, "verticalAerosolsIce")
       endif
!       if (.not. pkgTracerZSalinityActive) then                            !echmod deprecate
!          call finalize_stand_in_tracer_array(block, "verticalSalinity")
!       endif

       ! snow density tracer
       if (.not. pkgColumnTracerEffectiveSnowDensityActive) then
          call finalize_stand_in_tracer_array(block, "snowDensity")
       endif

       ! snow grain radius
       if (.not. pkgColumnTracerSnowGrainRadiusActive) then
          call finalize_stand_in_tracer_array(block, "snowGrainRadius")
          call finalize_stand_in_tracer_array(block, "snowLiquidMass")
          call finalize_stand_in_tracer_array(block, "snowIceMass")
       endif
       !-----------------------------------------------------------------------
       ! other column packages
       !-----------------------------------------------------------------------

       ! form drag
       call MPAS_pool_get_package(block % packages, "pkgColumnFormDragActive", pkgColumnFormDragActive)

       if (.not. pkgColumnFormDragActive) then

          call MPAS_pool_get_subpool(block % structs, "drag", drag)

          call MPAS_pool_get_field(drag, "oceanDragCoefficientSkin", oceanDragCoefficientSkin)
          call MPAS_pool_get_field(drag, "oceanDragCoefficientFloe", oceanDragCoefficientFloe)
          call MPAS_pool_get_field(drag, "oceanDragCoefficientKeel", oceanDragCoefficientKeel)
          call MPAS_pool_get_field(drag, "airDragCoefficientSkin", airDragCoefficientSkin)
          call MPAS_pool_get_field(drag, "airDragCoefficientFloe", airDragCoefficientFloe)
          call MPAS_pool_get_field(drag, "airDragCoefficientPond", airDragCoefficientPond)
          call MPAS_pool_get_field(drag, "airDragCoefficientRidge", airDragCoefficientRidge)
          call MPAS_pool_get_field(drag, "dragFreeboard", dragFreeboard)
          call MPAS_pool_get_field(drag, "dragIceSnowDraft", dragIceSnowDraft)
          call MPAS_pool_get_field(drag, "dragRidgeHeight", dragRidgeHeight)
          call MPAS_pool_get_field(drag, "dragRidgeSeparation", dragRidgeSeparation)
          call MPAS_pool_get_field(drag, "dragKeelDepth", dragKeelDepth)
          call MPAS_pool_get_field(drag, "dragKeelSeparation", dragKeelSeparation)
          call MPAS_pool_get_field(drag, "dragFloeLength", dragFloeLength)
          call MPAS_pool_get_field(drag, "dragFloeSeparation", dragFloeSeparation)

          oceanDragCoefficientSkin % array     => null()
          oceanDragCoefficientFloe % array     => null()
          oceanDragCoefficientKeel % array     => null()
          airDragCoefficientSkin % array       => null()
          airDragCoefficientFloe % array       => null()
          airDragCoefficientPond % array       => null()
          airDragCoefficientRidge % array      => null()
          dragFreeboard % array                => null()
          dragIceSnowDraft % array             => null()
          dragRidgeHeight % array              => null()
          dragRidgeSeparation % array          => null()
          dragKeelDepth % array                => null()
          dragKeelSeparation % array           => null()
          dragFloeLength % array               => null()
          dragFloeSeparation % array           => null()

       endif

       block => block % next
    end do

  end subroutine finalize_icepack_non_activated_pointers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  set_stand_in_tracer_array
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 5th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine set_stand_in_tracer_array(block, tracerName)

    type(block_type) :: block

    character(len=*), intent(in) :: &
         tracerName

    type(MPAS_pool_type), pointer :: &
         tracers

    type(field3DReal), pointer :: &
         tracerArray, &
         iceAreaCategory

    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_field(tracers, trim(tracerName), tracerArray, 1)
    call MPAS_pool_get_field(tracers, "iceAreaCategory", iceAreaCategory, 1)

    tracerArray % array => iceAreaCategory % array

  end subroutine set_stand_in_tracer_array

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  finalize_stand_in_tracer_array
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 29th October 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine finalize_stand_in_tracer_array(block, tracerName)

    type(block_type) :: block

    character(len=*), intent(in) :: &
         tracerName

    type(MPAS_pool_type), pointer :: &
         tracers

    type(field3DReal), pointer :: &
         tracerArray, &
         iceAreaCategory

    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_field(tracers, trim(tracerName), tracerArray, 1)

    tracerArray % array => null()

  end subroutine finalize_stand_in_tracer_array

!-----------------------------------------------------------------------
! other initialization
!-----------------------------------------------------------------------
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_history_variables
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_history_variables(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         ridging

    real(kind=RKIND), dimension(:,:), pointer :: &
         ratioRidgeThicknessToIce

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "ridging", ridging)
       call MPAS_pool_get_array(ridging, "ratioRidgeThicknessToIce", ratioRidgeThicknessToIce)

       ratioRidgeThicknessToIce = 1.0_RKIND

       block => block % next
    end do

  end subroutine init_column_history_variables

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack_initial_air_drag_coefficient
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  function seaice_icepack_initial_air_drag_coefficient() result(airDragCoefficient)

    use seaice_constants, only: &
         seaiceVonKarmanConstant, &
         seaiceIceSurfaceRoughness, &
         seaiceStabilityReferenceHeight

    real(kind=RKIND) :: airDragCoefficient

    ! atmo drag for RASM
    airDragCoefficient = (seaiceVonKarmanConstant/log(seaiceStabilityReferenceHeight/seaiceIceSurfaceRoughness)) &
                       * (seaiceVonKarmanConstant/log(seaiceStabilityReferenceHeight/seaiceIceSurfaceRoughness))

  end function seaice_icepack_initial_air_drag_coefficient

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack_reinitialize_fluxes
!
!> \brief
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_icepack_reinitialize_fluxes(domain)

    type(domain_type) :: domain

    ! atmospheric fluxes
    call seaice_column_reinitialize_atmospheric_fluxes(domain)

    ! oceanic fluxes
    call seaice_column_reinitialize_oceanic_fluxes(domain)

  end subroutine seaice_icepack_reinitialize_fluxes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_column_reinitialize_atmospheric_fluxes
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 31st August 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_column_reinitialize_atmospheric_fluxes(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         atmosFluxesPool, &
         shortwavePool, &
         atmosCouplingPool

    real(kind=RKIND), dimension(:), pointer :: &
         airStressCellU, &
         airStressCellV, &
         sensibleHeatFlux, &
         latentHeatFlux, &
         evaporativeWaterFlux, &
         longwaveUp, &
         absorbedShortwaveFlux, &
         atmosReferenceTemperature2m, &
         atmosReferenceHumidity2m, &
         atmosReferenceSpeed10m

    logical, pointer :: &
         config_use_column_physics

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCouplingPool)

       call MPAS_pool_get_array(velocitySolverPool, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressCellV", airStressCellV)

       call MPAS_pool_get_array(atmosCouplingPool, "atmosReferenceTemperature2m", atmosReferenceTemperature2m)
       call MPAS_pool_get_array(atmosCouplingPool, "atmosReferenceHumidity2m", atmosReferenceHumidity2m)
       call MPAS_pool_get_array(atmosCouplingPool, "atmosReferenceSpeed10m", atmosReferenceSpeed10m)

       airStressCellU(:)              = 0.0_RKIND
       airStressCellV(:)              = 0.0_RKIND

       atmosReferenceTemperature2m(:) = 0.0_RKIND
       atmosReferenceHumidity2m(:)    = 0.0_RKIND
       atmosReferenceSpeed10m(:)      = 0.0_RKIND

       call MPAS_pool_get_config(block % configs, "config_use_column_physics", config_use_column_physics)

       if (config_use_column_physics) then

          call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxesPool)
          call MPAS_pool_get_subpool(block % structs, "shortwave", shortwavePool)

          call MPAS_pool_get_array(atmosFluxesPool, "sensibleHeatFlux", sensibleHeatFlux)
          call MPAS_pool_get_array(atmosFluxesPool, "latentHeatFlux", latentHeatFlux)
          call MPAS_pool_get_array(atmosFluxesPool, "evaporativeWaterFlux", evaporativeWaterFlux)
          call MPAS_pool_get_array(atmosFluxesPool, "longwaveUp", longwaveUp)

          call MPAS_pool_get_array(shortwavePool, "absorbedShortwaveFlux", absorbedShortwaveFlux)

          absorbedShortwaveFlux(:) = 0.0_RKIND

          sensibleHeatFlux(:)      = 0.0_RKIND
          latentHeatFlux(:)        = 0.0_RKIND
          evaporativeWaterFlux(:)  = 0.0_RKIND
          longwaveUp(:)            = 0.0_RKIND

       endif ! config_use_column_physics

       block => block % next
    end do

  end subroutine seaice_column_reinitialize_atmospheric_fluxes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_column_reinitialize_oceanic_fluxes
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 31st August 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_column_reinitialize_oceanic_fluxes(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         oceanFluxesPool, &
         snowPool

    real(kind=RKIND), dimension(:), pointer :: &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         oceanShortwaveFlux, &
         snowLossToLeads, &
         snowMeltMassCell

    logical, pointer :: &
         config_use_column_physics, &
         config_use_column_snow_tracers

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_column_physics", config_use_column_physics)
       call MPAS_pool_get_config(block % configs, "config_use_column_snow_tracers", config_use_column_snow_tracers)

       if (config_use_column_physics) then

          call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxesPool)

          call MPAS_pool_get_array(oceanFluxesPool, "oceanFreshWaterFlux", oceanFreshWaterFlux)
          call MPAS_pool_get_array(oceanFluxesPool, "oceanSaltFlux", oceanSaltFlux)
          call MPAS_pool_get_array(oceanFluxesPool, "oceanHeatFlux", oceanHeatFlux)
          call MPAS_pool_get_array(oceanFluxesPool, "oceanShortwaveFlux", oceanShortwaveFlux)

          oceanFreshWaterFlux(:) = 0.0_RKIND
          oceanSaltFlux(:)       = 0.0_RKIND
          oceanHeatFlux(:)       = 0.0_RKIND
          oceanShortwaveFlux(:)  = 0.0_RKIND

          if (config_use_column_snow_tracers) then
             call MPAS_pool_get_subpool(block % structs, "snow", snowPool)

             call MPAS_pool_get_array(snowPool, "snowLossToLeads", snowLossToLeads)
             call MPAS_pool_get_array(snowPool, "snowMeltMassCell", snowMeltMassCell)

             snowLossToLeads(:) = 0.0_RKIND
             snowMeltMassCell(:) = 0.0_RKIND

          endif ! config_use_column_snow_tracers

       endif ! config_use_column_physics

       block => block % next
    end do

  end subroutine seaice_column_reinitialize_oceanic_fluxes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_tracer_object_for_biogeochemistry
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 14 September 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_tracer_object_for_biogeochemistry(domain, tracerObject)

    use icepack_intfc, only: icepack_init_zbgc
    use icepack_intfc, only: icepack_init_parameters
    use icepack_intfc, only: icepack_query_parameters

    type(domain_type), intent(in) :: &
         domain

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    logical, pointer :: &
         config_use_brine, &
!         config_use_vertical_zsalinity, &      !echmod deprecate
         config_use_vertical_biochemistry, &
         config_use_vertical_tracers, &
         config_use_skeletal_biochemistry, &
         config_use_shortwave_bioabsorption, &
         config_use_nitrate, &
         config_use_carbon, &
         config_use_chlorophyll, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_nonreactive, &
         config_use_humics, &
         config_use_DON, &
         config_use_iron, &
         config_use_zaerosols

    real(kind=RKIND), pointer :: &
         config_new_ice_fraction_biotracer, &
         config_fraction_biotracer_in_frazil, &
         config_ratio_Si_to_N_diatoms, &
         config_ratio_Si_to_N_small_plankton, &
         config_ratio_Si_to_N_phaeocystis, &
         config_ratio_S_to_N_diatoms, &
         config_ratio_S_to_N_small_plankton, &
         config_ratio_S_to_N_phaeocystis, &
         config_ratio_Fe_to_C_diatoms, &
         config_ratio_Fe_to_C_small_plankton, &
         config_ratio_Fe_to_C_phaeocystis, &
         config_ratio_Fe_to_N_diatoms, &
         config_ratio_Fe_to_N_small_plankton, &
         config_ratio_Fe_to_N_phaeocystis, &
         config_ratio_Fe_to_DON, &
         config_ratio_Fe_to_DOC_saccharids, &
         config_ratio_Fe_to_DOC_lipids, &
         config_chla_absorptivity_of_diatoms, &
         config_chla_absorptivity_of_small_plankton, &
         config_chla_absorptivity_of_phaeocystis, &
         config_light_attenuation_diatoms, &
         config_light_attenuation_small_plankton, &
         config_light_attenuation_phaeocystis, &
         config_light_inhibition_diatoms, &
         config_light_inhibition_small_plankton, &
         config_light_inhibition_phaeocystis, &
         config_maximum_growth_rate_diatoms, &
         config_maximum_growth_rate_small_plankton, &
         config_maximum_growth_rate_phaeocystis, &
         config_temperature_growth_diatoms, &
         config_temperature_growth_small_plankton, &
         config_temperature_growth_phaeocystis, &
         config_grazed_fraction_diatoms, &
         config_grazed_fraction_small_plankton, &
         config_grazed_fraction_phaeocystis, &
         config_mortality_diatoms, &
         config_mortality_small_plankton, &
         config_mortality_phaeocystis, &
         config_temperature_mortality_diatoms, &
         config_temperature_mortality_small_plankton, &
         config_temperature_mortality_phaeocystis, &
         config_exudation_diatoms, &
         config_exudation_small_plankton, &
         config_exudation_phaeocystis, &
         config_nitrate_saturation_diatoms, &
         config_nitrate_saturation_small_plankton, &
         config_nitrate_saturation_phaeocystis, &
         config_ammonium_saturation_diatoms, &
         config_ammonium_saturation_small_plankton, &
         config_ammonium_saturation_phaeocystis, &
         config_silicate_saturation_diatoms, &
         config_silicate_saturation_small_plankton, &
         config_silicate_saturation_phaeocystis, &
         config_iron_saturation_diatoms, &
         config_iron_saturation_small_plankton, &
         config_iron_saturation_phaeocystis, &
         config_fraction_spilled_to_DON, &
         config_degredation_of_DON, &
         config_fraction_DON_ammonium, &
         config_fraction_loss_to_saccharids, &
         config_fraction_loss_to_lipids, &
         config_fraction_exudation_to_saccharids, &
         config_fraction_exudation_to_lipids, &
         config_remineralization_saccharids, &
         config_remineralization_lipids, &
         config_mobility_type_diatoms, &
         config_mobility_type_small_plankton, &
         config_mobility_type_phaeocystis, &
         config_mobility_type_saccharids, &
         config_mobility_type_lipids, &
         config_mobility_type_inorganic_carbon, &
         config_mobility_type_proteins, &
         config_mobility_type_dissolved_iron, &
         config_mobility_type_particulate_iron, &
         config_mobility_type_black_carbon1, &
         config_mobility_type_black_carbon2, &
         config_mobility_type_dust1, &
         config_mobility_type_dust2, &
         config_mobility_type_dust3, &
         config_mobility_type_dust4, &
         config_ratio_C_to_N_diatoms, &
         config_ratio_C_to_N_small_plankton, &
         config_ratio_C_to_N_phaeocystis, &
         config_ratio_chla_to_N_diatoms, &
         config_ratio_chla_to_N_small_plankton, &
         config_ratio_chla_to_N_phaeocystis, &
         config_scales_absorption_diatoms, &
         config_scales_absorption_small_plankton, &
         config_scales_absorption_phaeocystis, &
         config_ratio_C_to_N_proteins, &
         config_mobility_type_nitrate, &
         config_mobility_type_ammonium, &
         config_mobility_type_DMSPp, &
         config_mobility_type_DMSPd, &
         config_mobility_type_silicate, &
         config_mobility_type_humics, &
         config_rapid_mobile_to_stationary_time, &
         config_long_mobile_to_stationary_time

    integer, pointer :: &
         ONE, &
         nIceLayers, &
         nSnowLayers, &
         nBioLayers, &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON, &
         nParticulateIron, &
         nDissolvedIron, &
         nzAerosols, &
         nZBGCTracers, &
         maxAerosolType, &
         maxAlgaeType, &
         maxDOCType, &
         maxDICType, &
         maxDONType, &
         maxIronType

    logical :: &
         use_nitrogen

    integer :: &
         nTracers_temp, &
         iAerosols

    real(kind=RKIND) :: &
         rtmp1, &
         rtmp2
    real (kind=RKIND), dimension (:), allocatable :: &
        BGCTracerType, &
        initialMobileFraction, &
        retentionTime, &
        releaseTime, &
        newIceBGCFraction

    ! save tracer array size
    nTracers_temp = tracerObject % nTracers
    tracerObject % nTracers = tracerObject % nTracersNotBio

    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)
!    call MPAS_pool_get_config(domain % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)  !echmod deprecate
    call MPAS_pool_get_config(domain % configs, "config_use_shortwave_bioabsorption", config_use_shortwave_bioabsorption)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_tracers", config_use_vertical_tracers)
    call MPAS_pool_get_config(domain % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_nitrate", config_use_nitrate)
    call MPAS_pool_get_config(domain % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(domain % configs, "config_use_chlorophyll", config_use_chlorophyll)
    call MPAS_pool_get_config(domain % configs, "config_use_ammonium", config_use_ammonium)
    call MPAS_pool_get_config(domain % configs, "config_use_silicate", config_use_silicate)
    call MPAS_pool_get_config(domain % configs, "config_use_DMS", config_use_DMS)
    call MPAS_pool_get_config(domain % configs, "config_use_nonreactive", config_use_nonreactive)
    call MPAS_pool_get_config(domain % configs, "config_use_humics", config_use_humics)
    call MPAS_pool_get_config(domain % configs, "config_use_DON", config_use_DON)
    call MPAS_pool_get_config(domain % configs, "config_use_iron", config_use_iron)
    call MPAS_pool_get_config(domain % configs, "config_use_zaerosols", config_use_zaerosols)
    call MPAS_pool_get_config(domain % configs, "config_new_ice_fraction_biotracer", config_new_ice_fraction_biotracer)
    call MPAS_pool_get_config(domain % configs, "config_fraction_biotracer_in_frazil", config_fraction_biotracer_in_frazil)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Si_to_N_diatoms", config_ratio_Si_to_N_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Si_to_N_small_plankton", config_ratio_Si_to_N_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Si_to_N_phaeocystis", config_ratio_Si_to_N_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_ratio_S_to_N_diatoms", config_ratio_S_to_N_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_ratio_S_to_N_small_plankton", config_ratio_S_to_N_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_ratio_S_to_N_phaeocystis", config_ratio_S_to_N_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_C_diatoms", config_ratio_Fe_to_C_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_C_small_plankton", config_ratio_Fe_to_C_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_C_phaeocystis", config_ratio_Fe_to_C_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_N_diatoms", config_ratio_Fe_to_N_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_N_small_plankton", config_ratio_Fe_to_N_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_N_phaeocystis", config_ratio_Fe_to_N_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_DON", config_ratio_Fe_to_DON)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_DOC_saccharids", config_ratio_Fe_to_DOC_saccharids)
    call MPAS_pool_get_config(domain % configs, "config_ratio_Fe_to_DOC_lipids", config_ratio_Fe_to_DOC_lipids)
    call MPAS_pool_get_config(domain % configs, "config_rapid_mobile_to_stationary_time", config_rapid_mobile_to_stationary_time)
    call MPAS_pool_get_config(domain % configs, "config_long_mobile_to_stationary_time", config_long_mobile_to_stationary_time)
    call MPAS_pool_get_config(domain % configs, "config_chla_absorptivity_of_diatoms", config_chla_absorptivity_of_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_chla_absorptivity_of_small_plankton", &
                                                 config_chla_absorptivity_of_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_chla_absorptivity_of_phaeocystis", config_chla_absorptivity_of_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_light_attenuation_diatoms", config_light_attenuation_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_light_attenuation_small_plankton", config_light_attenuation_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_light_attenuation_phaeocystis", config_light_attenuation_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_light_inhibition_diatoms", config_light_inhibition_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_light_inhibition_small_plankton", config_light_inhibition_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_light_inhibition_phaeocystis", config_light_inhibition_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_maximum_growth_rate_diatoms", config_maximum_growth_rate_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_maximum_growth_rate_small_plankton", &
                                                 config_maximum_growth_rate_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_maximum_growth_rate_phaeocystis", config_maximum_growth_rate_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_temperature_growth_diatoms", config_temperature_growth_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_temperature_growth_small_plankton", &
                                                 config_temperature_growth_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_temperature_growth_phaeocystis", config_temperature_growth_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_grazed_fraction_diatoms", config_grazed_fraction_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_grazed_fraction_small_plankton", config_grazed_fraction_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_grazed_fraction_phaeocystis", config_grazed_fraction_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_mortality_diatoms", config_mortality_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_mortality_small_plankton", config_mortality_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_mortality_phaeocystis", config_mortality_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_temperature_mortality_diatoms", config_temperature_mortality_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_temperature_mortality_small_plankton", &
                                                 config_temperature_mortality_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_temperature_mortality_phaeocystis", &
                                                 config_temperature_mortality_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_exudation_diatoms", config_exudation_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_exudation_small_plankton", config_exudation_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_exudation_phaeocystis", config_exudation_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_nitrate_saturation_diatoms", config_nitrate_saturation_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_nitrate_saturation_small_plankton", &
                                                 config_nitrate_saturation_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_nitrate_saturation_phaeocystis", config_nitrate_saturation_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_ammonium_saturation_diatoms", config_ammonium_saturation_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_ammonium_saturation_small_plankton", &
                                                 config_ammonium_saturation_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_ammonium_saturation_phaeocystis", config_ammonium_saturation_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_silicate_saturation_diatoms", config_silicate_saturation_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_silicate_saturation_small_plankton", &
                                                 config_silicate_saturation_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_silicate_saturation_phaeocystis", config_silicate_saturation_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_iron_saturation_diatoms", config_iron_saturation_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_iron_saturation_small_plankton", config_iron_saturation_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_iron_saturation_phaeocystis", config_iron_saturation_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_fraction_spilled_to_DON", config_fraction_spilled_to_DON)
    call MPAS_pool_get_config(domain % configs, "config_degredation_of_DON", config_degredation_of_DON)
    call MPAS_pool_get_config(domain % configs, "config_fraction_DON_ammonium", config_fraction_DON_ammonium)
    call MPAS_pool_get_config(domain % configs, "config_fraction_loss_to_saccharids", config_fraction_loss_to_saccharids)
    call MPAS_pool_get_config(domain % configs, "config_fraction_loss_to_lipids",  config_fraction_loss_to_lipids)
    call MPAS_pool_get_config(domain % configs, "config_fraction_exudation_to_saccharids", config_fraction_exudation_to_saccharids)
    call MPAS_pool_get_config(domain % configs, "config_fraction_exudation_to_lipids", config_fraction_exudation_to_lipids)
    call MPAS_pool_get_config(domain % configs, "config_remineralization_saccharids", config_remineralization_saccharids)
    call MPAS_pool_get_config(domain % configs, "config_remineralization_lipids", config_remineralization_lipids)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_diatoms", config_mobility_type_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_small_plankton", config_mobility_type_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_phaeocystis", config_mobility_type_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_nitrate", config_mobility_type_nitrate)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_ammonium", config_mobility_type_ammonium)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_silicate", config_mobility_type_silicate)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_DMSPp", config_mobility_type_DMSPp)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_DMSPd", config_mobility_type_DMSPd)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_humics", config_mobility_type_humics)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_saccharids", config_mobility_type_saccharids)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_lipids", config_mobility_type_lipids)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_inorganic_carbon", config_mobility_type_inorganic_carbon)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_proteins", config_mobility_type_proteins)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_dissolved_iron", config_mobility_type_dissolved_iron)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_particulate_iron", config_mobility_type_particulate_iron)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_black_carbon1", config_mobility_type_black_carbon1)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_black_carbon2", config_mobility_type_black_carbon2)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_dust1", config_mobility_type_dust1)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_dust2", config_mobility_type_dust2)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_dust3", config_mobility_type_dust3)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_dust4", config_mobility_type_dust4)
    call MPAS_pool_get_config(domain % configs, "config_ratio_C_to_N_diatoms", config_ratio_C_to_N_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_ratio_C_to_N_small_plankton", config_ratio_C_to_N_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_ratio_C_to_N_phaeocystis", config_ratio_C_to_N_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_ratio_chla_to_N_diatoms", config_ratio_chla_to_N_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_ratio_chla_to_N_small_plankton", config_ratio_chla_to_N_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_ratio_chla_to_N_phaeocystis", config_ratio_chla_to_N_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_scales_absorption_diatoms", config_scales_absorption_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_scales_absorption_small_plankton", config_scales_absorption_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_scales_absorption_phaeocystis", config_scales_absorption_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_ratio_C_to_N_proteins", config_ratio_C_to_N_proteins)

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "ONE", ONE)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nSnowLayers", nSnowLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nBioLayers",nBioLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nAlgae", nAlgae)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDOC", nDOC)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDIC", nDIC)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDON", nDON)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nParticulateIron", nParticulateIron)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDissolvedIron", nDissolvedIron)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nzAerosols", nzAerosols)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nZBGCTracers", nZBGCTracers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxAerosolType", maxAerosolType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxAlgaeType", maxAlgaeType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxDOCType", maxDOCType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxDICType", maxDICType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxDONType", maxDONType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxIronType", maxIronType)

    allocate(BGCTracerType(nZBGCTracers))
    allocate(initialMobileFraction(nZBGCTracers))
    allocate(retentionTime(nZBGCTracers))
    allocate(releaseTime(nZBGCTracers))
    allocate(newIceBGCFraction(nZBGCTracers))

    retentionTime(:) = 0.0_RKIND
    releaseTime(:) = 0.0_RKIND
    initialMobileFraction(:) = 0.0_RKIND
    BGCTracerType(:) = 0.0_RKIND
    newIceBGCFraction(:) = 0.0_RKIND

    use_nitrogen = .false.
    if (config_use_skeletal_biochemistry .or. config_use_vertical_biochemistry) &
        use_nitrogen = .true.

    allocate(tracerObject % index_verticalAerosolsConc(maxAerosolType))
    allocate(tracerObject % index_verticalAerosolsConcLayer(maxAerosolType))
    allocate(tracerObject % index_verticalAerosolsConcShortwave(maxAerosolType))
    tracerObject % nzAerosolsIndex = nzAerosols

    allocate(tracerObject % index_algaeConc(maxAlgaeType))
    allocate(tracerObject % index_algaeConcLayer(maxAlgaeType))
    tracerObject % nAlgaeIndex = nAlgae

    allocate(tracerObject % index_algalCarbon(maxAlgaeType))
    allocate(tracerObject % index_algalCarbonLayer(maxAlgaeType))
    tracerObject % nAlgalCarbonIndex = nAlgae

    allocate(tracerObject % index_DOCConc(maxDOCType))
    allocate(tracerObject % index_DOCConcLayer(maxDOCType))
    tracerObject % nDOCIndex = nDOC

    allocate(tracerObject % index_DICConc(maxDICType))
    allocate(tracerObject % index_DICConcLayer(maxDICType))
    tracerObject % nDICIndex = nDIC

    allocate(tracerObject % index_algalChlorophyll(maxAlgaeType))
    allocate(tracerObject % index_algalChlorophyllLayer(maxAlgaeType))
    tracerObject % nAlgalChlorophyllIndex = nAlgae

    allocate(tracerObject % index_DONConc(maxDONType))
    allocate(tracerObject % index_DONConcLayer(maxDONType))
    tracerObject % nDONIndex = nDON

    allocate(tracerObject % index_particulateIronConc(maxIronType))
    allocate(tracerObject % index_particulateIronConcLayer(maxIronType))
    tracerObject % nParticulateIronIndex = nParticulateIron

    allocate(tracerObject % index_dissolvedIronConc(maxIronType))
    allocate(tracerObject % index_dissolvedIronConcLayer(maxIronType))
    tracerObject % nDissolvedIronIndex = nDissolvedIron

    call icepack_init_parameters(&
         ratio_Si2N_diatoms_in=config_ratio_Si_to_N_diatoms, &
         ratio_Si2N_sp_in=config_ratio_Si_to_N_small_plankton, &
         ratio_Si2N_phaeo_in=config_ratio_Si_to_N_phaeocystis, &
         ratio_S2N_diatoms_in=config_ratio_S_to_N_diatoms, &
         ratio_S2N_sp_in=config_ratio_S_to_N_small_plankton, &
         ratio_S2N_phaeo_in=config_ratio_S_to_N_phaeocystis, &
         ratio_Fe2C_diatoms_in=config_ratio_Fe_to_C_diatoms, &
         ratio_Fe2C_sp_in=config_ratio_Fe_to_C_small_plankton, &
         ratio_Fe2C_phaeo_in=config_ratio_Fe_to_C_phaeocystis, &
         ratio_Fe2N_diatoms_in=config_ratio_Fe_to_N_diatoms, &
         ratio_Fe2N_sp_in=config_ratio_Fe_to_N_small_plankton, &
         ratio_Fe2N_phaeo_in=config_ratio_Fe_to_N_phaeocystis, &
         ratio_Fe2DON_in=config_ratio_Fe_to_DON, &
         ratio_Fe2DOC_s_in=config_ratio_Fe_to_DOC_saccharids, &
         ratio_Fe2DOC_l_in=config_ratio_Fe_to_DOC_lipids, &
         chlabs_diatoms_in=config_chla_absorptivity_of_diatoms, &
         chlabs_sp_in=config_chla_absorptivity_of_small_plankton, &
         chlabs_phaeo_in=config_chla_absorptivity_of_phaeocystis, &
         alpha2max_low_diatoms_in=config_light_attenuation_diatoms, &
         alpha2max_low_sp_in=config_light_attenuation_small_plankton, &
         alpha2max_low_phaeo_in=config_light_attenuation_phaeocystis, &
         beta2max_diatoms_in=config_light_inhibition_diatoms, &
         beta2max_sp_in=config_light_inhibition_small_plankton, &
         beta2max_phaeo_in=config_light_inhibition_phaeocystis, &
         mu_max_diatoms_in=config_maximum_growth_rate_diatoms, &
         mu_max_sp_in=config_maximum_growth_rate_small_plankton, &
         mu_max_phaeo_in=config_maximum_growth_rate_phaeocystis, &
         grow_Tdep_diatoms_in=config_temperature_growth_diatoms, &
         grow_Tdep_sp_in=config_temperature_growth_small_plankton, &
         grow_Tdep_phaeo_in=config_temperature_growth_phaeocystis, &
         fr_graze_diatoms_in=config_grazed_fraction_diatoms, &
         fr_graze_sp_in=config_grazed_fraction_small_plankton, &
         fr_graze_phaeo_in=config_grazed_fraction_phaeocystis, &
         mort_pre_diatoms_in=config_mortality_diatoms, &
         mort_pre_sp_in=config_mortality_small_plankton, &
         mort_pre_phaeo_in=config_mortality_phaeocystis, &
         mort_Tdep_diatoms_in=config_temperature_mortality_diatoms, &
         mort_Tdep_sp_in=config_temperature_mortality_small_plankton, &
         mort_Tdep_phaeo_in=config_temperature_mortality_phaeocystis, &
         k_exude_diatoms_in=config_exudation_diatoms, &
         k_exude_sp_in=config_exudation_small_plankton, &
         k_exude_phaeo_in=config_exudation_phaeocystis, &
         K_Nit_diatoms_in=config_nitrate_saturation_diatoms, &
         K_Nit_sp_in=config_nitrate_saturation_small_plankton, &
         K_Nit_phaeo_in=config_nitrate_saturation_phaeocystis, &
         K_Am_diatoms_in=config_ammonium_saturation_diatoms, &
         K_Am_sp_in=config_ammonium_saturation_small_plankton, &
         K_Am_phaeo_in=config_ammonium_saturation_phaeocystis, &
         K_Sil_diatoms_in=config_silicate_saturation_diatoms, &
         K_Sil_sp_in=config_silicate_saturation_small_plankton, &
         K_Sil_phaeo_in=config_silicate_saturation_phaeocystis, &
         K_Fe_diatoms_in=config_iron_saturation_diatoms, &
         K_Fe_sp_in=config_iron_saturation_small_plankton, &
         K_Fe_phaeo_in=config_iron_saturation_phaeocystis, &
         f_don_protein_in=config_fraction_spilled_to_DON, &
         kn_bac_protein_in=config_degredation_of_DON, &
         f_don_Am_protein_in=config_fraction_DON_ammonium, &
         f_doc_s_in=config_fraction_loss_to_saccharids, &
         f_doc_l_in=config_fraction_loss_to_lipids, &
         f_exude_s_in=config_fraction_exudation_to_saccharids, &
         f_exude_l_in=config_fraction_exudation_to_lipids, &
         k_bac_s_in=config_remineralization_saccharids, &
         k_bac_l_in=config_remineralization_lipids, &
         algaltype_diatoms_in=config_mobility_type_diatoms, &
         algaltype_sp_in=config_mobility_type_small_plankton, &
         algaltype_phaeo_in=config_mobility_type_phaeocystis, &
         doctype_s_in=config_mobility_type_saccharids, &
         doctype_l_in=config_mobility_type_lipids, &
         dictype_1_in=config_mobility_type_inorganic_carbon, &
         dontype_protein_in=config_mobility_type_proteins, &
         fedtype_1_in=config_mobility_type_dissolved_iron, &
         feptype_1_in=config_mobility_type_particulate_iron, &
         zaerotype_bc1_in=config_mobility_type_black_carbon1, &
         zaerotype_bc2_in=config_mobility_type_black_carbon2, &
         zaerotype_dust1_in=config_mobility_type_dust1, &
         zaerotype_dust2_in=config_mobility_type_dust2, &
         zaerotype_dust3_in=config_mobility_type_dust3, &
         zaerotype_dust4_in=config_mobility_type_dust4, &
         ratio_C2N_diatoms_in=config_ratio_C_to_N_diatoms, &
         ratio_C2N_sp_in=config_ratio_C_to_N_small_plankton, &
         ratio_C2N_phaeo_in=config_ratio_C_to_N_phaeocystis, &
         ratio_chl2N_diatoms_in=config_ratio_chla_to_N_diatoms, &
         ratio_chl2N_sp_in=config_ratio_chla_to_N_small_plankton, &
         ratio_chl2N_phaeo_in=config_ratio_chla_to_N_phaeocystis, &
         F_abs_chl_diatoms_in=config_scales_absorption_diatoms, &
         F_abs_chl_sp_in=config_scales_absorption_small_plankton, &
         F_abs_chl_phaeo_in=config_scales_absorption_phaeocystis, &
         ratio_C2N_proteins_in=config_ratio_C_to_N_proteins, &
         nitratetype_in=config_mobility_type_nitrate, &
         ammoniumtype_in=config_mobility_type_ammonium, &
         dmspptype_in=config_mobility_type_DMSPp, &
         dmspdtype_in=config_mobility_type_DMSPd, &
         silicatetype_in=config_mobility_type_silicate, &
         humtype_in=config_mobility_type_humics, &
         frazil_scav_in=config_fraction_biotracer_in_frazil, &
         initbio_frac_in=config_new_ice_fraction_biotracer, &
         tau_min_in=config_rapid_mobile_to_stationary_time, &
         tau_max_in=config_long_mobile_to_stationary_time)

    call init_zbgc_tracer_indices(domain, tracerObject, use_nitrogen, nTracers_temp, &
                                 BGCTracerType, initialMobileFraction, retentionTime, &
                                 releaseTime, newIceBGCFraction)

    call icepack_init_zbgc(bgc_tracer_type_in=BGCTracerType,        &
                           zbgc_frac_init_in=initialMobileFraction, &
                           tau_ret_in=retentionTime,                &
                           tau_rel_in=releaseTime,                  &
                           zbgc_init_frac_in=newIceBGCFraction)

    ! check calculated tracer array size
    if (nTracers_temp /= tracerObject % nTracers) then
       call mpas_log_write(&
            "init_column_tracer_object_for_biogeochemistry: nTracers_temp: $i, nTracers: $i", &
            messageType=MPAS_LOG_CRIT, intArgs=(/nTracers_temp, tracerObject % nTracers/))
    endif

    call mpas_log_write(" ----- compare values after init_zbgc -----")

    rtmp1 = config_new_ice_fraction_biotracer
    call icepack_query_parameters(initbio_frac_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('initbio_frac differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = config_rapid_mobile_to_stationary_time
    call icepack_query_parameters(tau_min_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('tau_min differs $r $r',realArgs=(/rtmp1,rtmp2/))

    rtmp1 = config_long_mobile_to_stationary_time
    call icepack_query_parameters(tau_max_out=rtmp2)
    if (rtmp1 /= rtmp2) call mpas_log_write('tau_max differs $r $r',realArgs=(/rtmp1,rtmp2/))

    deallocate(BGCTracerType)
    deallocate(initialMobileFraction)
    deallocate(retentionTime)
    deallocate(releaseTime)
    deallocate(newIceBGCFraction)

  end subroutine init_column_tracer_object_for_biogeochemistry

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_zbgc_tracer_indices
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 24 May 2024
!> \details
!> Moves icepack_init_zbgc_tracer_indices to the interface
!
!-----------------------------------------------------------------------

  subroutine init_zbgc_tracer_indices(domain, tracerObject, use_nitrogen, nTracersTemp, &
                                      BGCTracerType, initialMobileFraction, retentionTime, &
                                      releaseTime, newIceBGCFraction)

    type(domain_type), intent(in) :: &
         domain

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    real(kind=RKIND), dimension(:), intent(inout) :: &
        BGCTracerType, &
        initialMobileFraction, &
        retentionTime, &
        releaseTime, &
        newIceBGCFraction

    integer, intent(in) :: &
         nTracersTemp

    logical, intent(in) :: &
         use_nitrogen

    logical, pointer :: &
         config_use_brine, &
!         config_use_vertical_zsalinity, &      !echmod deprecate
         config_use_vertical_biochemistry, &
         config_use_vertical_tracers, &
         config_use_skeletal_biochemistry, &
         config_use_shortwave_bioabsorption, &
         config_use_nitrate, &
         config_use_carbon, &
         config_use_chlorophyll, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_nonreactive, &
         config_use_humics, &
         config_use_DON, &
         config_use_iron, &
         config_use_zaerosols

    real(kind=RKIND), pointer :: &
         config_mobility_type_diatoms, &
         config_mobility_type_small_plankton, &
         config_mobility_type_phaeocystis, &
         config_mobility_type_nitrate, &
         config_mobility_type_ammonium, &
         config_mobility_type_silicate, &
         config_mobility_type_DMSPp, &
         config_mobility_type_DMSPd, &
         config_mobility_type_humics, &
         config_mobility_type_saccharids, &
         config_mobility_type_lipids, &
         config_mobility_type_inorganic_carbon, &
         config_mobility_type_proteins, &
         config_mobility_type_dissolved_iron, &
         config_mobility_type_particulate_iron, &
         config_mobility_type_black_carbon1, &
         config_mobility_type_black_carbon2, &
         config_mobility_type_dust1, &
         config_mobility_type_dust2, &
         config_mobility_type_dust3, &
         config_mobility_type_dust4, &
         config_rapid_mobile_to_stationary_time, &
         config_long_mobile_to_stationary_time, &
         config_fraction_biotracer_in_frazil, &
         config_new_ice_fraction_biotracer

    integer, pointer :: &
         nIceLayers, &
         nSnowLayers, &
         nBioLayers, &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON, &
         nParticulateIron, &
         nDissolvedIron, &
         nzAerosols, &
         nZBGCTracers, &
         maxAerosolType, &
         maxAlgaeType, &
         maxDOCType, &
         maxDICType, &
         maxDONType, &
         maxIronType

    integer :: &
         iAerosols, &
         nTracerDependOption, &
         nCount, &
         nTracerDepend, &
         iBioLayer, &
         iBioTracer

    real(kind=RKIND) :: &
         rtmp1, &
         rtmp2

    real (kind=RKIND), dimension (:), allocatable :: &
        algalType, &
        docType, &
        dicType, &
        donType, &
        fedType, &
        fepType, &
        zAeroType

    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)
!    call MPAS_pool_get_config(domain % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)  !echmod deprecate
    call MPAS_pool_get_config(domain % configs, "config_use_shortwave_bioabsorption", config_use_shortwave_bioabsorption)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_tracers", config_use_vertical_tracers)
    call MPAS_pool_get_config(domain % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_nitrate", config_use_nitrate)
    call MPAS_pool_get_config(domain % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(domain % configs, "config_use_chlorophyll", config_use_chlorophyll)
    call MPAS_pool_get_config(domain % configs, "config_use_ammonium", config_use_ammonium)
    call MPAS_pool_get_config(domain % configs, "config_use_silicate", config_use_silicate)
    call MPAS_pool_get_config(domain % configs, "config_use_DMS", config_use_DMS)
    call MPAS_pool_get_config(domain % configs, "config_use_nonreactive", config_use_nonreactive)
    call MPAS_pool_get_config(domain % configs, "config_use_humics", config_use_humics)
    call MPAS_pool_get_config(domain % configs, "config_use_DON", config_use_DON)
    call MPAS_pool_get_config(domain % configs, "config_use_iron", config_use_iron)
    call MPAS_pool_get_config(domain % configs, "config_use_zaerosols", config_use_zaerosols)

    call MPAS_pool_get_config(domain % configs, "config_mobility_type_diatoms", config_mobility_type_diatoms)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_small_plankton", config_mobility_type_small_plankton)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_phaeocystis", config_mobility_type_phaeocystis)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_nitrate", config_mobility_type_nitrate)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_ammonium", config_mobility_type_ammonium)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_silicate", config_mobility_type_silicate)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_DMSPp", config_mobility_type_DMSPp)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_DMSPd", config_mobility_type_DMSPd)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_humics", config_mobility_type_humics)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_saccharids", config_mobility_type_saccharids)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_lipids", config_mobility_type_lipids)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_inorganic_carbon", config_mobility_type_inorganic_carbon)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_proteins", config_mobility_type_proteins)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_dissolved_iron", config_mobility_type_dissolved_iron)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_particulate_iron", config_mobility_type_particulate_iron)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_black_carbon1", config_mobility_type_black_carbon1)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_black_carbon2", config_mobility_type_black_carbon2)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_dust1", config_mobility_type_dust1)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_dust2", config_mobility_type_dust2)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_dust3", config_mobility_type_dust3)
    call MPAS_pool_get_config(domain % configs, "config_mobility_type_dust4", config_mobility_type_dust4)
    call MPAS_pool_get_config(domain % configs, "config_rapid_mobile_to_stationary_time", &
                                                 config_rapid_mobile_to_stationary_time)
    call MPAS_pool_get_config(domain % configs, "config_long_mobile_to_stationary_time", &
                                                 config_long_mobile_to_stationary_time)
    call MPAS_pool_get_config(domain % configs, "config_fraction_biotracer_in_frazil", &
                                                 config_fraction_biotracer_in_frazil)
    call MPAS_pool_get_config(domain % configs, "config_new_ice_fraction_biotracer", &
                                                 config_new_ice_fraction_biotracer)

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nSnowLayers", nSnowLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nBioLayers",nBioLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nAlgae", nAlgae)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDOC", nDOC)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDIC", nDIC)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDON", nDON)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nParticulateIron", nParticulateIron)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nDissolvedIron", nDissolvedIron)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nzAerosols", nzAerosols)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nZBGCTracers", nZBGCTracers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxAerosolType", maxAerosolType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxAlgaeType", maxAlgaeType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxDOCType", maxDOCType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxDICType", maxDICType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxDONType", maxDONType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxIronType", maxIronType)

      tracerObject % nTracersNotBio = tracerObject % nTracers
      tracerObject % index_brineFraction = 0
      if (config_use_brine) then
          tracerObject % index_brineFraction = tracerObject % nTracers + 1   ! ice volume fraction with salt
          tracerObject % nTracers = tracerObject % nTracers + 1
          tracerObject % parentIndex(tracerObject % index_brineFraction)   = 1   ! volume-weighted
          tracerObject % firstAncestorMask  (tracerObject % index_brineFraction,1) = 0.0_RKIND  ! volume-weighted
          tracerObject % firstAncestorMask  (tracerObject % index_brineFraction,2) = 1.0_RKIND  ! volume-weighted
          tracerObject % firstAncestorMask  (tracerObject % index_brineFraction,3) = 0.0_RKIND  ! volume-weighted
          tracerObject % ancestorNumber(tracerObject % index_brineFraction) = 0
          tracerObject % ancestorIndices  (tracerObject % index_brineFraction,1) = 0
          tracerObject % ancestorIndices  (tracerObject % index_brineFraction,2) = 0
      endif

      nTracerDependOption = 0                    ! if tracerObject % index_brineFraction /= 0 then use fbri dependency
      if (tracerObject % index_brineFraction == 0) nTracerDependOption = -1 ! otherwise make tracers depend on ice volume

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

      tracerObject % nBioTracers = 0
      tracerObject % nBioTracersShortwave = 0

      ! vectors of size maxAlgaeType
      tracerObject % index_algaeConcLayer(:) = 0
      tracerObject % index_algalCarbonLayer(:) = 0
      tracerObject % index_algalChlorophyllLayer(:) = 0
      tracerObject % index_algaeConc(:) = 0
      tracerObject % index_algalCarbon(:) = 0
      tracerObject % index_algalChlorophyll(:) = 0

      ! vectors of size maxDICType
      tracerObject % index_DICConcLayer(:) = 0
      tracerObject % index_DICConc(:) = 0

      ! vectors of size maxDOCType
      tracerObject % index_DOCConcLayer(:) = 0
      tracerObject % index_DOCConc(:) = 0

      ! vectors of size maxDONType
      tracerObject % index_DONConcLayer(:) = 0
      tracerObject % index_DONConc(:) = 0

      ! vectors of size maxIronType
      tracerObject % index_dissolvedIronConcLayer(:) = 0
      tracerObject % index_particulateIronConcLayer(:) = 0
      tracerObject % index_dissolvedIronConc(:) = 0
      tracerObject % index_particulateIronConc(:) = 0

      ! vectors of size maxAerosolType
      tracerObject % index_verticalAerosolsConcLayer(:) = 0
      tracerObject % index_verticalAerosolsConcShortwave(:) = 0
      tracerObject % index_verticalAerosolsConc(:) = 0

      tracerObject % index_nitrateConcLayer    = 0
      tracerObject % index_ammoniumConcLayer   = 0
      tracerObject % index_silicateConcLayer    = 0
      tracerObject % index_DMSPpConcLayer  = 0
      tracerObject % index_DMSPdConcLayer  = 0
      tracerObject % index_DMSConcLayer    = 0
      tracerObject % index_nonreactiveConcLayer = 0
      tracerObject % index_humicsConcLayer = 0
      tracerObject % index_chlorophyllShortwave     = 0
      tracerObject % index_LayerIndexToBioIndex(:)  = 0
      tracerObject % index_LayerIndexToDataArray(:) = 0

      tracerObject % index_nitrateConc    = 0
      tracerObject % index_ammoniumConc   = 0
      tracerObject % index_silicateConc   = 0
      tracerObject % index_DMSPpConc  = 0
      tracerObject % index_DMSPdConc  = 0
      tracerObject % index_DMSConc    = 0
      tracerObject % index_nonreactiveConc    = 0
      tracerObject % index_humicsConc    = 0

      allocate(algalType(maxAlgaeType))
      allocate(docType(maxDOCType))
      allocate(dicType(maxDICType))
      allocate(donType(maxDONType))
      allocate(fedType(maxIronType))
      allocate(fepType(maxIronType))
      allocate(zAeroType(maxAerosolType))

      algalType(1) = config_mobility_type_diatoms
      algalType(2) = config_mobility_type_small_plankton
      algalType(3) = config_mobility_type_phaeocystis

      docType(1) = config_mobility_type_saccharids
      docType(2) = config_mobility_type_lipids

      dicType(1) = config_mobility_type_inorganic_carbon

      donType(1) = config_mobility_type_proteins

      fedType(1) = config_mobility_type_dissolved_iron
      fepType(1) = config_mobility_type_particulate_iron

      zAeroType(1) = config_mobility_type_black_carbon1
      zAeroType(2) = config_mobility_type_black_carbon2
      zAeroType(3) = config_mobility_type_dust1
      zAeroType(4) = config_mobility_type_dust2
      zAeroType(5) = config_mobility_type_dust3
      zAeroType(6) = config_mobility_type_dust4

      if (config_use_skeletal_biochemistry) then

         nCount = 1
         nTracerDepend = 0

         if (config_use_shortwave_bioabsorption) then
           tracerObject % index_chlorophyllShortwave = 1
           tracerObject % nBioTracersShortwave = nIceLayers+nSnowLayers+2  ! only the bottom layer
                                                                           ! will be nonzero
          endif
      elseif (config_use_vertical_tracers) then ! defined on nBioLayers+1 in ice
                                                ! and 2 snow layers (snow surface + interior)
         nCount = nBioLayers + 1
         nTracerDepend = 2 + tracerObject % index_brineFraction + nTracerDependOption

         if (use_nitrogen) then
            if (config_use_shortwave_bioabsorption) then
               tracerObject % index_chlorophyllShortwave = 1
               tracerObject % nBioTracersShortwave = nIceLayers+nSnowLayers+2
            endif
         endif ! use_nitrogen

      endif ! config_use_skeletal_biochemistry or config_use_vertical_tracers

      if (config_use_skeletal_biochemistry .or. config_use_vertical_tracers) then

      !-----------------------------------------------------------------
      ! assign tracer indices and dependencies
      ! BGCTracerType: < 0  purely mobile , >= 0 stationary
      !------------------------------------------------------------------

      if (use_nitrogen) then
         do iBioTracer = 1, nAlgae
            call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction,  &
                                      tracerObject % index_algaeConc(iBioTracer),     &
                                      tracerObject % index_algaeConcLayer(iBioTracer),&
                                      algalType(iBioTracer),  nTracerDepend,          &
                                      tracerObject % nTracers,                        &
                                      tracerObject % nBioTracers,                     &
                                      BGCTracerType,                                  &
                                      tracerObject % parentIndex,                     &
                                      tracerObject % firstAncestorMask,               &
                                      tracerObject % ancestorNumber,                  &
                                      tracerObject % ancestorIndices,                 &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_algaeConcLayer(iBioTracer)) &
                           = iBioTracer
         enddo   ! iBioTracer
      endif ! use_nitrogen

      if (config_use_nitrate) then
         call init_bgc_tracer_indices(&
                                      nCount, tracerObject % index_brineFraction,  &
                                      tracerObject % index_nitrateConc,            &
                                      tracerObject % index_nitrateConcLayer,       &
                                      config_mobility_type_nitrate, nTracerDepend, &
                                      tracerObject % nTracers,                     &
                                      tracerObject % nBioTracers,                  &
                                      BGCTracerType,                               &
                                      tracerObject % parentIndex,                  &
                                      tracerObject % firstAncestorMask,            &
                                      tracerObject % ancestorNumber,               &
                                      tracerObject % ancestorIndices,              &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_nitrateConcLayer) &
                           = maxAlgaeType + 1
      endif ! config_use_nitrate

      if (config_use_carbon) then
       !
       ! Algal C is not yet distinct from algal N
       ! * Reqires exudation and/or changing C:N ratios
       ! for implementation
       !
       !  do iBioTracer = 1,nAlgae
       !     call init_bgc_tracer_indices(nCount,tracerObject % index_brineFraction,      &
       !                               tracerObject % index_algalCarbon(iBioTracer),      &
       !                               tracerObject % index_algalCarbonLayer(iBioTracer), &
       !                               algalType(iBioTracer),   nTracerDepend,            &
       !                               tracerObject % nTracers,                           &
       !                               tracerObject % nBioTracers,                        &
       !                               BGCTracerType,                                     &
       !                               tracerObject % parentIndex,                        &
       !                               tracerObject % firstAncestorMask,                  &
       !                               tracerObject % ancestorNumber,                     &
       !                               tracerObject % ancestorIndices,                    &
       !                               tracerObject % index_LayerIndexToBioIndex)
       !     tracerObject % index_LayerIndexToDataArray(tracerObject % index_algalCarbonLayer(iBioTracer)) &
       !                    = maxAlgaeType + 1 + iBioTracer
       !  enddo   ! iBioTracer

         do iBioTracer = 1, nDOC
            call init_bgc_tracer_indices(&
                                      nCount, tracerObject % index_brineFraction,    &
                                      tracerObject % index_DOCConc(iBioTracer),      &
                                      tracerObject % index_DOCConcLayer(iBioTracer), &
                                      docType(iBioTracer), nTracerDepend,            &
                                      tracerObject % nTracers,                       &
                                      tracerObject % nBioTracers,                    &
                                      BGCTracerType,                                 &
                                      tracerObject % parentIndex,                    &
                                      tracerObject % firstAncestorMask,              &
                                      tracerObject % ancestorNumber,                 &
                                      tracerObject % ancestorIndices,                &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_DOCConcLayer(iBioTracer)) &
                           = maxAlgaeType + 1 + iBioTracer
         enddo   ! iBioTracer
         do iBioTracer = 1, nDIC
            call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction,  &
                                      tracerObject % index_DICConc(iBioTracer),       &
                                      tracerObject % index_DICConcLayer(iBioTracer),  &
                                      dicType(iBioTracer), nTracerDepend,             &
                                      tracerObject % nTracers,                        &
                                      tracerObject % nBioTracers,                     &
                                      BGCTracerType,                                  &
                                      tracerObject % parentIndex,                     &
                                      tracerObject % firstAncestorMask,               &
                                      tracerObject % ancestorNumber,                  &
                                      tracerObject % ancestorIndices,                 &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_DICConcLayer(iBioTracer)) &
                           = maxAlgaeType + maxDOCType + 1 + iBioTracer
         enddo   ! iBioTracer
      endif      ! config_use_carbon

      if (config_use_chlorophyll) then
         do iBioTracer = 1, nAlgae
            call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction,          &
                                      tracerObject % index_algalChlorophyll(iBioTracer),      &
                                      tracerObject % index_algalChlorophyllLayer(iBioTracer), &
                                      algalType(iBioTracer), nTracerDepend,                   &
                                      tracerObject % nTracers,                                &
                                      tracerObject % nBioTracers,                             &
                                      BGCTracerType,                                          &
                                      tracerObject % parentIndex,                             &
                                      tracerObject % firstAncestorMask,                       &
                                      tracerObject % ancestorNumber,                          &
                                      tracerObject % ancestorIndices,                         &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_algalChlorophyllLayer(iBioTracer)) &
                           = maxAlgaeType + 1 + maxDOCType + maxDICType + iBioTracer
         enddo   ! iBioTracer
      endif      ! config_use_chlorophyll

      if (config_use_ammonium) then
         call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction, &
                                      tracerObject % index_ammoniumConc,          &
                                      tracerObject % index_ammoniumConcLayer,     &
                                      config_mobility_type_ammonium,              &
                                      nTracerDepend,                              &
                                      tracerObject % nTracers,                    &
                                      tracerObject % nBioTracers,                 &
                                      BGCTracerType,                              &
                                      tracerObject % parentIndex,                 &
                                      tracerObject % firstAncestorMask,           &
                                      tracerObject % ancestorNumber,              &
                                      tracerObject % ancestorIndices,             &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_ammoniumConcLayer) &
                           = 2*maxAlgaeType + maxDOCType + maxDICType + 2
      endif
      if (config_use_silicate) then
         call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction, &
                                      tracerObject % index_silicateConc,          &
                                      tracerObject % index_silicateConcLayer,     &
                                      config_mobility_type_silicate,              &
                                      nTracerDepend,                              &
                                      tracerObject % nTracers,                    &
                                      tracerObject % nBioTracers,                 &
                                      BGCTracerType,                              &
                                      tracerObject % parentIndex,                 &
                                      tracerObject % firstAncestorMask,           &
                                      tracerObject % ancestorNumber,              &
                                      tracerObject % ancestorIndices,             &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_silicateConcLayer) &
                           = 2*maxAlgaeType + maxDOCType + maxDICType + 3
      endif
      if (config_use_DMS) then   ! all together
         call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction, &
                                      tracerObject % index_DMSPpConc,             &
                                      tracerObject % index_DMSPpConcLayer,        &
                                      config_mobility_type_DMSPp, nTracerDepend,  &
                                      tracerObject % nTracers,                    &
                                      tracerObject % nBioTracers,                 &
                                      BGCTracerType,                              &
                                      tracerObject % parentIndex,                 &
                                      tracerObject % firstAncestorMask,           &
                                      tracerObject % ancestorNumber,              &
                                      tracerObject % ancestorIndices,             &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_DMSPpConcLayer) &
                           = 2*maxAlgaeType + maxDOCType + maxDICType + 4

            call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction, &
                                      tracerObject % index_DMSPdConc,                &
                                      tracerObject % index_DMSPdConcLayer,           &
                                      config_mobility_type_DMSPd, nTracerDepend,     &
                                      tracerObject % nTracers,                       &
                                      tracerObject % nBioTracers,                    &
                                      BGCTracerType,                                 &
                                      tracerObject % parentIndex,                    &
                                      tracerObject % firstAncestorMask,              &
                                      tracerObject % ancestorNumber,                 &
                                      tracerObject % ancestorIndices,                &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_DMSPdConcLayer) &
                           = 2*maxAlgaeType + maxDOCType + maxDICType + 5

            call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction, &
                                      tracerObject % index_DMSConc,                  &
                                      tracerObject % index_DMSConcLayer,             &
                                      config_mobility_type_DMSPd, nTracerDepend,     &
                                      tracerObject % nTracers,                       &
                                      tracerObject % nBioTracers,                    &
                                      BGCTracerType,                                 &
                                      tracerObject % parentIndex,                    &
                                      tracerObject % firstAncestorMask,              &
                                      tracerObject % ancestorNumber,                 &
                                      tracerObject % ancestorIndices,                &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_DMSConcLayer) &
                           = 2*maxAlgaeType + maxDOCType + maxDICType + 6
      endif
      if (config_use_nonreactive) then
         call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction, &
                                      tracerObject % index_nonreactiveConc,       &
                                      tracerObject % index_nonreactiveConcLayer,  &
                                      config_mobility_type_nitrate, nTracerDepend,&
                                      tracerObject % nTracers,                    &
                                      tracerObject % nBioTracers,                 &
                                      BGCTracerType,                              &
                                      tracerObject % parentIndex,                 &
                                      tracerObject % firstAncestorMask,           &
                                      tracerObject % ancestorNumber,              &
                                      tracerObject % ancestorIndices,             &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_nonreactiveConcLayer) &
                           = 2*maxAlgaeType + maxDOCType + maxDICType + 7
      endif
      if (config_use_DON) then
         do iBioTracer = 1, nDON
            call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction, &
                                      tracerObject % index_DONConc(iBioTracer),      &
                                      tracerObject % index_DONConcLayer(iBioTracer), &
                                      donType(iBioTracer), nTracerDepend,            &
                                      tracerObject % nTracers,                       &
                                      tracerObject % nBioTracers,                    &
                                      BGCTracerType,                                 &
                                      tracerObject % parentIndex,                    &
                                      tracerObject % firstAncestorMask,              &
                                      tracerObject % ancestorNumber,                 &
                                      tracerObject % ancestorIndices,                &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_DONConcLayer(iBioTracer)) &
                           = 2*maxAlgaeType + maxDOCType + maxDICType + 7 + iBioTracer
         enddo   ! iBioTracer
      endif      ! config_use_DON
      if (config_use_iron) then
         do iBioTracer = 1, nDissolvedIron
            call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction,          &
                                      tracerObject % index_dissolvedIronConc(iBioTracer),     &
                                      tracerObject % index_dissolvedIronConcLayer(iBioTracer),&
                                      fedType(iBioTracer), nTracerDepend,                     &
                                      tracerObject % nTracers,                                &
                                      tracerObject % nBioTracers,                             &
                                      BGCTracerType,   tracerObject % parentIndex,            &
                                      tracerObject % firstAncestorMask,                       &
                                      tracerObject % ancestorNumber,                          &
                                      tracerObject % ancestorIndices,                         &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_dissolvedIronConcLayer(iBioTracer)) &
                           = 2*maxAlgaeType + maxDOCType + maxDICType + maxDONType + 7 + iBioTracer
         enddo   ! iBioTracer
         do iBioTracer = 1, nParticulateIron
            call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction,             &
                                      tracerObject % index_particulateIronConc(iBioTracer),      &
                                      tracerObject % index_particulateIronConcLayer(iBioTracer), &
                                      fepType(iBioTracer), nTracerDepend,                        &
                                      tracerObject % nTracers,                                   &
                                      tracerObject % nBioTracers,                                &
                                      BGCTracerType,   tracerObject % parentIndex,               &
                                      tracerObject % firstAncestorMask,                          &
                                      tracerObject % ancestorNumber,                             &
                                      tracerObject % ancestorIndices,                            &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_particulateIronConcLayer(iBioTracer)) &
                           = 2*maxAlgaeType + maxDOCType + maxDICType + maxDONType + maxIronType + 7 + iBioTracer
         enddo   ! iBioTracer
      endif      ! config_use_iron

      if (config_use_humics) then
         call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction, &
                                      tracerObject % index_humicsConc,            &
                                      tracerObject % index_humicsConcLayer,       &
                                      config_mobility_type_humics, nTracerDepend, &
                                      tracerObject % nTracers,                    &
                                      tracerObject % nBioTracers,                 &
                                      BGCTracerType,                              &
                                      tracerObject % parentIndex,                 &
                                      tracerObject % firstAncestorMask,           &
                                      tracerObject % ancestorNumber,              &
                                      tracerObject % ancestorIndices,             &
                                      tracerObject % index_LayerIndexToBioIndex)
            tracerObject % index_LayerIndexToDataArray(tracerObject % index_humicsConcLayer) &
                           = 2*maxAlgaeType + maxDOCType + 8 + maxDICType + maxDONType + 2*maxIronType + maxAerosolType
      endif
      endif  ! config_use_skeletal_biochemistry or config_use_vertical_tracers

      if (config_use_vertical_tracers) then ! defined on nBioLayers+1 in ice
                                            ! and 2 snow layers (snow surface + interior)
         nCount = nBioLayers + 1
         nTracerDepend = 2 + tracerObject % index_brineFraction + nTracerDependOption

         ! z layer aerosols
         if (config_use_zaerosols) then
            do iBioTracer = 1, nzAerosols
               if (config_use_shortwave_bioabsorption) then
                  tracerObject % index_verticalAerosolsConcShortwave(iBioTracer) = &
                                 tracerObject % nBioTracersShortwave + 1
                  tracerObject % nBioTracersShortwave = &
                                 tracerObject % nBioTracersShortwave + nIceLayers + nSnowLayers+2
               endif
               call init_bgc_tracer_indices(nCount, tracerObject % index_brineFraction,              &
                                         tracerObject % index_verticalAerosolsConc(iBioTracer),      &
                                         tracerObject % index_verticalAerosolsConcLayer(iBioTracer), &
                                         zAeroType(iBioTracer), nTracerDepend,                       &
                                         tracerObject % nTracers,                                    &
                                         tracerObject % nBioTracers,                                 &
                                         BGCTracerType, tracerObject % parentIndex,                  &
                                         tracerObject % firstAncestorMask,                           &
                                         tracerObject % ancestorNumber,                              &
                                         tracerObject % ancestorIndices,                             &
                                         tracerObject % index_LayerIndexToBioIndex)
               tracerObject % index_LayerIndexToDataArray(tracerObject % index_verticalAerosolsConcLayer(iBioTracer)) &
                              = 2*maxAlgaeType + maxDOCType + maxDICType + maxDONType + 2*maxIronType + 7 + iBioTracer
            enddo   ! iBioTracer
         endif      ! config_use_zaerosols

         tracerObject % index_mobileFraction = 0
         if (tracerObject % nBioTracers > 0) then
            tracerObject % index_mobileFraction = tracerObject % nTracers + 1
            tracerObject % nTracers = tracerObject % nTracers + tracerObject % nBioTracers
            do iBioLayer = 1,tracerObject % nBioTracers
               initialMobileFraction(iBioLayer) = 1.0_RKIND
               tracerObject % parentIndex(tracerObject % index_mobileFraction+iBioLayer-1) &
                              = 2+tracerObject % index_brineFraction
               tracerObject % firstAncestorMask(tracerObject % index_mobileFraction+ iBioLayer - 1,1)  = 0.0_RKIND
               tracerObject % firstAncestorMask(tracerObject % index_mobileFraction+ iBioLayer - 1,2)  = 1.0_RKIND
               tracerObject % firstAncestorMask(tracerObject % index_mobileFraction+ iBioLayer - 1,3)  = 0.0_RKIND
               tracerObject % ancestorNumber(tracerObject % index_mobileFraction+ iBioLayer - 1) = 1
               tracerObject % ancestorIndices(tracerObject % index_mobileFraction+ iBioLayer - 1,1) &
                              = tracerObject % index_brineFraction
               tracerObject % ancestorIndices(tracerObject % index_mobileFraction+ iBioLayer - 1,2) = 0
               retentionTime(iBioLayer) = 1.0_RKIND
               releaseTime(iBioLayer) = 1.0_RKIND

               if (BGCTracerType(iBioLayer) >=  0.0_RKIND .and. BGCTracerType(iBioLayer) < 0.5_RKIND) then
                  retentionTime(iBioLayer) = config_rapid_mobile_to_stationary_time
                  releaseTime(iBioLayer) = config_long_mobile_to_stationary_time
                  initialMobileFraction(iBioLayer) = 1.0_RKIND
               elseif (BGCTracerType(iBioLayer) >= 0.5_RKIND .and. BGCTracerType(iBioLayer) < 1.0_RKIND) then
                  retentionTime(iBioLayer) = config_rapid_mobile_to_stationary_time
                  releaseTime(iBioLayer) = config_rapid_mobile_to_stationary_time
                  initialMobileFraction(iBioLayer) = 1.0_RKIND
               elseif (BGCTracerType(iBioLayer) >= 1.0_RKIND .and. BGCTracerType(iBioLayer) < 2.0_RKIND) then
                  retentionTime(iBioLayer) = config_long_mobile_to_stationary_time
                  releaseTime(iBioLayer) = config_rapid_mobile_to_stationary_time
                  initialMobileFraction(iBioLayer) = 1.0_RKIND
               elseif (BGCTracerType(iBioLayer) >= 2.0_RKIND ) then
                  retentionTime(iBioLayer) = config_long_mobile_to_stationary_time
                  releaseTime(iBioLayer) = config_long_mobile_to_stationary_time
                  initialMobileFraction(iBioLayer) = 1.0_RKIND
               endif
            enddo
         endif

      endif ! config_use_vertical_tracers

      do iBioLayer = 1, tracerObject % nBioTracers
         newIceBGCFraction(iBioLayer) = config_fraction_biotracer_in_frazil
         if (BGCTracerType(iBioLayer) < 0.0_RKIND) &
               newIceBGCFraction(iBioLayer) = config_new_ice_fraction_biotracer
      enddo

      if (.NOT.config_use_shortwave_bioabsorption) tracerObject % nBioTracersShortwave = 1

      !-----------------------------------------------------------------
      ! final consistency checks
      !-----------------------------------------------------------------

      if (tracerObject % nBioTracers > nZBGCTracers) then
         call mpas_log_write(&
            "init_zbgc_tracer_indices: nBioTracers $i, greater than max nZBGCTracers: $i", &
            messageType=MPAS_LOG_CRIT, intArgs=(/tracerObject % nBioTracers, nZBGCTracers/))
      endif

      deallocate(algalType)
      deallocate(docType)
      deallocate(dicType)
      deallocate(donType)
      deallocate(fedType)
      deallocate(fepType)
      deallocate(zAeroType)

  end subroutine init_zbgc_tracer_indices

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_bgc_tracer_indices
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 24 May 2024
!> \details
!>
!> Moves icepack_init_bgc_trcr to the interface
!
!-----------------------------------------------------------------------

  subroutine init_bgc_tracer_indices( nCount,          iBrineFraction, &
                                      iBGCConc,        iBGCConcLayer,  &
                                      BGCType,         nTracerDepend,  &
                                      nTracers,        nBioTracers,    &
                                      BGCTracerType,   parentIndexBGC, &
                                      firstAncestorMaskBGC,            &
                                      ancestorNumberBGC,               &
                                      ancestorIndicesBGC,              &
                                      layerIndexToBioIndex)

      integer, intent(in) :: &
         nCount           , & ! counter
         nTracerDepend    , & ! tracer dependency index
         iBrineFraction

      integer, intent(inout) :: &
         nTracers         , & ! number of tracers
         nBioTracers      , & ! number of bio tracers
         iBGCConc         , & ! tracer index
         iBGCConcLayer        ! bio tracer index

      integer, dimension(:), intent(inout) :: &
         parentIndexBGC   , & ! tracer dependencies
         ancestorNumberBGC, & ! number of underlying tracer layers
         layerIndexToBioIndex

      integer, dimension(:,:), intent(inout) :: &
         ancestorIndicesBGC    ! indices of underlying tracer layers

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         firstAncestorMaskBGC  ! = 0 or 1 depending on tracer dependency
                               ! argument 2:  (1) aice, (2) vice, (3) vsno

      real (kind=RKIND), intent(in) :: &
         BGCType               ! bio tracer transport type (mobile vs stationary)

      real (kind=RKIND), dimension(:), intent(inout) :: &
         BGCTracerType         ! bio tracer transport type array

      ! local variables

      integer :: &
         iCount , & ! loop index
         nTemp  , & ! temporary values
         ancestorIndicesBGC1, & ! temporary values
         ancestorIndicesBGC2

      real (kind=RKIND) :: &
         firstAncestorMaskBGC1, & ! temporary values
         firstAncestorMaskBGC2, &
         firstAncestorMaskBGC3

         iBGCConc = nTracers + 1
         nBioTracers = nBioTracers + 1
         iBGCConcLayer = nBioTracers
         BGCTracerType(nBioTracers) = BGCType

         if (nCount > 1) then
            ! include vertical bgc in snow
            do iCount = nCount, nCount+1
               nTracers = nTracers + 1
               parentIndexBGC  (iBGCConc + iCount  ) = 2 ! snow volume
               firstAncestorMaskBGC    (iBGCConc + iCount,1) = 0.0_RKIND
               firstAncestorMaskBGC    (iBGCConc + iCount,2) = 0.0_RKIND
               firstAncestorMaskBGC    (iBGCConc + iCount,3) = 1.0_RKIND
               ancestorNumberBGC(iBGCConc + iCount  ) = 0
               ancestorIndicesBGC    (iBGCConc + iCount,1) = 0
               ancestorIndicesBGC    (iBGCConc + iCount,2) = 0
            enddo

            firstAncestorMaskBGC1 = 0.0_RKIND
            firstAncestorMaskBGC2 = 1.0_RKIND
            firstAncestorMaskBGC3 = 0.0_RKIND
            nTemp = 1
            ancestorIndicesBGC1 = iBrineFraction
            ancestorIndicesBGC2 = 0
         else  ! nCount = 1
            firstAncestorMaskBGC1 = 1.0_RKIND
            firstAncestorMaskBGC2 = 0.0_RKIND
            firstAncestorMaskBGC3 = 0.0_RKIND
            nTemp = 0
            ancestorIndicesBGC1 = 0
            ancestorIndicesBGC2 = 0
         endif ! nCount

         do iCount = 1, nCount     !in ice
            nTracers = nTracers + 1
            parentIndexBGC  (iBGCConc + iCount - 1  ) = nTracerDepend
            firstAncestorMaskBGC    (iBGCConc + iCount - 1,1) = firstAncestorMaskBGC1
            firstAncestorMaskBGC    (iBGCConc + iCount - 1,2) = firstAncestorMaskBGC2
            firstAncestorMaskBGC    (iBGCConc + iCount - 1,3) = firstAncestorMaskBGC3
            ancestorNumberBGC(iBGCConc + iCount - 1  ) = nTemp
            ancestorIndicesBGC    (iBGCConc + iCount - 1,1) = ancestorIndicesBGC1
            ancestorIndicesBGC    (iBGCConc + iCount - 1,2) = ancestorIndicesBGC2
         enddo

         layerIndexToBioIndex (iBGCConcLayer) = iBGCConc

  end subroutine init_bgc_tracer_indices

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_biogeochemistry
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 17th September 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_biogeochemistry_profiles(domain, tracerObject)

    use icepack_intfc, only: &
         icepack_init_bgc, &
         icepack_init_hbrine, &
         icepack_load_ocean_bio_array

    type(domain_type), intent(inout) :: domain

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         biogeochemistry, &
         diagnostics_biogeochemistry, &
         ocean_coupling, &
         tracers

    logical, pointer :: &
         config_use_brine, &
         config_use_vertical_tracers, &
         config_use_skeletal_biochemistry, &
         config_do_restart_hbrine

    real(kind=RKIND), pointer :: &
         config_dt, &
         config_snow_porosity_at_ice_surface

    real(kind=RKIND), dimension(:), pointer :: &
         oceanNitrateConc, &
         oceanSilicateConc, &
         oceanAmmoniumConc, &
         oceanDMSConc, &
         oceanDMSPConc, &
         oceanHumicsConc, &
         seaSurfaceSalinity , &
         verticalGrid, &          ! cgrid
         interfaceBiologyGrid, &  ! igrid
         biologyGrid, &           ! bgrid
         verticalShortwaveGrid, & ! swgrid
         interfaceGrid, &         ! icgrid
         DOCPoolFractions

    real(kind=RKIND), dimension(:,:), pointer :: &
         oceanAlgaeConc, &
         oceanDOCConc, &
         oceanDICConc, &
         oceanDONConc, &
         oceanParticulateIronConc, &
         oceanDissolvedIronConc, &
         oceanZAerosolConc, &
         oceanBioConcentrations, &
         totalVerticalBiologyIce, &
         totalVerticalBiologySnow

    integer, dimension(:,:), pointer :: &
         newlyFormedIce

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         bioPorosity, &
         bioDiffusivity, &
         bioTemperature, &
         bioPermeability, &
         bioShortwaveFlux, &
         bioTracerShortwave, &
         iceSalinity, &
         brineFraction

    integer, pointer :: &
         nCellsSolve

    integer :: &
         iCell

    logical :: &
         abortFlag, &
         setGetPhysicsTracers, &
         setGetBGCTracers

    character(len=strKIND) :: &
         abortMessage

    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)
    call MPAS_pool_get_config(domain % configs, "config_do_restart_hbrine", config_do_restart_hbrine)
    call MPAS_pool_get_config(domain % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_tracers", config_use_vertical_tracers)
    call MPAS_pool_get_config(domain % configs, "config_dt", config_dt)
    call MPAS_pool_get_config(domain % configs, "config_snow_porosity_at_ice_surface", config_snow_porosity_at_ice_surface)

    abortFlag = .false.

    setGetPhysicsTracers = .false.
    setGetBGCTracers     = .true.

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)
       call MPAS_pool_get_subpool(block % structs, "diagnostics_biogeochemistry", diagnostics_biogeochemistry)

       call MPAS_pool_get_array(biogeochemistry, "bioPorosity", bioPorosity)
       call MPAS_pool_get_array(biogeochemistry, "bioDiffusivity", bioDiffusivity)
       call MPAS_pool_get_array(biogeochemistry, "bioTemperature", bioTemperature)
       call MPAS_pool_get_array(biogeochemistry, "bioPermeability", bioPermeability)
       call MPAS_pool_get_array(biogeochemistry, "bioShortwaveFlux", bioShortwaveFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioConcentrations", oceanBioConcentrations)
       call MPAS_pool_get_array(biogeochemistry, "totalVerticalBiologyIce", totalVerticalBiologyIce)
       call MPAS_pool_get_array(biogeochemistry, "totalVerticalBiologySnow", totalVerticalBiologySnow)

       call MPAS_pool_get_array(biogeochemistry, "bioTracerShortwave", bioTracerShortwave)
       call MPAS_pool_get_array(biogeochemistry, "interfaceBiologyGrid", interfaceBiologyGrid)
       call MPAS_pool_get_array(biogeochemistry, "interfaceGrid", interfaceGrid)
       call MPAS_pool_get_array(biogeochemistry, "verticalGrid", verticalGrid)
       call MPAS_pool_get_array(biogeochemistry, "biologyGrid", biologyGrid)
       call MPAS_pool_get_array(biogeochemistry, "verticalShortwaveGrid", verticalShortwaveGrid)
       call MPAS_pool_get_array(biogeochemistry, "oceanAlgaeConc", oceanAlgaeConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanDOCConc", oceanDOCConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanDICConc", oceanDICConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanDONConc", oceanDONConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanParticulateIronConc", oceanParticulateIronConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanDissolvedIronConc", oceanDissolvedIronConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanNitrateConc", oceanNitrateConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanSilicateConc", oceanSilicateConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanAmmoniumConc", oceanAmmoniumConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSConc", oceanDMSConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSPConc", oceanDMSPConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanHumicsConc", oceanHumicsConc)
       call MPAS_pool_get_array(biogeochemistry, "oceanZAerosolConc", oceanZAerosolConc)
       call MPAS_pool_get_array(biogeochemistry, "newlyFormedIce", newlyFormedIce)
       call MPAS_pool_get_array(biogeochemistry, "DOCPoolFractions", DOCPoolFractions)

       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
       call MPAS_pool_get_array(ocean_coupling, "seaSurfaceSalinity", seaSurfaceSalinity)

       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity, 1)
       call MPAS_pool_get_array(tracers, "brineFraction", brineFraction, 1)

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       call icepack_init_hbrine(&
            bgrid_out=biologyGrid, &
            igrid_out=interfaceBiologyGrid, &
            cgrid_out=verticalGrid, &
            icgrid_out=interfaceGrid, &
            swgrid_out=verticalShortwaveGrid, &
            phi_snow=config_snow_porosity_at_ice_surface)

       do iCell = 1, nCellsSolve
          if (.not. config_do_restart_hbrine) then
             ! initialize newly formed ice
             newlyFormedIce(:,iCell) = 1

             ! initialize brine fraction
             brineFraction(:,:,iCell) = 1.0_RKIND
          endif

          ! set the category tracer array
          call set_cice_tracer_array_category(block, tracerObject, &
               tracerArrayCategory, iCell, setGetPhysicsTracers, setGetBGCTracers)

          if (config_use_vertical_tracers .or. config_use_skeletal_biochemistry) then

             call icepack_load_ocean_bio_array(&
                  nit=oceanNitrateConc(iCell), &
                  amm=oceanAmmoniumConc(iCell), &
                  sil=oceanSilicateConc(iCell), &
                  dmsp=oceanDMSPConc(iCell), &
                  dms=oceanDMSConc(iCell), &
                  algalN=oceanAlgaeConc(:,iCell), &
                  doc=oceanDOCConc(:,iCell), &
                  don=oceanDONConc(:,iCell), &
                  dic=oceanDICConc(:,iCell), &
                  fed=oceanDissolvedIronConc(:,iCell), &
                  fep=oceanParticulateIronConc(:,iCell), &
                  zaeros=oceanZAerosolConc(:,iCell), &
                  ocean_bio_all=oceanBioConcentrations(:,iCell), &
                  hum=oceanHumicsConc(iCell))

             call icepack_init_bgc(&
                  sicen=iceSalinity(:,:,iCell), &
                  trcrn=tracerArrayCategory(tracerObject % nTracersNotBio+1:tracerObject % nTracers,:), &
                  sss=seaSurfaceSalinity(iCell), &
                  ocean_bio_all=oceanBioConcentrations(:,iCell), &
                  DOCPoolFractions=DOCPoolFractions)

             call seaice_icepack_write_warnings(icepack_warnings_aborted())

          endif ! biogeochemistry

          ! get the category tracer array
          call get_cice_tracer_array_category(block, tracerObject, &
               tracerArrayCategory, iCell, setGetPhysicsTracers, setGetBGCTracers)

       enddo ! iCell

       block => block % next
    end do

  end subroutine init_column_biogeochemistry_profiles

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack_reinitialize_diagnostics_thermodynamics
!
!> \brief Reinitialize thermodynamics diagnostics
!> \author Adrian K. Turner, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Andrew Roberts, Erin Thomas, Jon Wolfe, LANL, Anthony Craig, NOAA (cntr), David Bailey, NCAR
!> \date 2022-2023
!> \details
!>  Reinitialize thermodynamics diagnostics
!
!-----------------------------------------------------------------------

  subroutine seaice_icepack_reinitialize_diagnostics_thermodynamics(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         atmosFluxesPool, &
         meltGrowthRatesPool, &
         diagnosticsPool, &
         tracersAggregatePool, &
         pondsPool, &
         shortwavePool, &
         dragPool, &
         snowPool

    ! atmospheric fluxes
    real(kind=RKIND), dimension(:), pointer :: &
         surfaceHeatFlux, &
         surfaceConductiveFlux

    real(kind=RKIND), dimension(:,:), pointer :: &
         surfaceHeatFluxCategory, &
         surfaceConductiveFluxCategory, &
         latentHeatFluxCategory, &
         sensibleHeatFluxCategory

    ! melt growth rates
    real(kind=RKIND), dimension(:), pointer :: &
         congelation, &
         frazilFormation, &
         snowiceFormation, &
         snowThicknessChange, &
         surfaceIceMelt, &
         snowMelt, &
         basalIceMelt, &
         lateralIceMelt

    ! snow model
    real(kind=RKIND), dimension(:), pointer :: &
         snowLossToLeads, &
         snowMeltMassCell, &
         snowDensityViaContent, &
         snowDensityViaCompaction, &
         snowRadiusInStandardRadiationScheme

    real(kind=RKIND), dimension(:,:), pointer :: &
         snowMeltMassCategory, &
         snowRadiusInStandardRadiationSchemeCategory

    ! diagnostic tendencies
    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaTendencyThermodynamics, &
         iceVolumeTendencyThermodynamics, &
         iceAgeTendencyThermodynamics, &
         iceAreaCell, &
         iceVolumeCell, &
         iceAgeCell

    ! pond fluxes
    real(kind=RKIND), dimension(:), pointer :: &
         pondFreshWaterFlux

    ! shortwave
    real(kind=RKIND), dimension(:), pointer :: &
         bareIceAlbedoCell, &
         snowAlbedoCell, &
         pondAlbedoCell

    ! drag variables
    real(kind=RKIND), pointer :: &
         config_ice_ocean_drag_coefficient

    real(kind=RKIND), dimension(:), pointer :: &
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
         dragFloeSeparation

    logical, pointer :: &
         config_use_ice_age, &
         config_use_form_drag, &
         config_use_column_physics, &
         config_use_column_snow_tracers

    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_physics", config_use_column_physics)

    if (config_use_column_physics) then

       block => domain % blocklist
       do while (associated(block))

          ! atmospheric fluxes
          call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxesPool)

          call MPAS_pool_get_array(atmosFluxesPool, "surfaceHeatFlux", surfaceHeatFlux)
          call MPAS_pool_get_array(atmosFluxesPool, "surfaceConductiveFlux", surfaceConductiveFlux)
          call MPAS_pool_get_array(atmosFluxesPool, "surfaceHeatFluxCategory", surfaceHeatFluxCategory)
          call MPAS_pool_get_array(atmosFluxesPool, "surfaceConductiveFluxCategory", surfaceConductiveFluxCategory)
          call MPAS_pool_get_array(atmosFluxesPool, "latentHeatFluxCategory", latentHeatFluxCategory)
          call MPAS_pool_get_array(atmosFluxesPool, "sensibleHeatFluxCategory", sensibleHeatFluxCategory)

          surfaceHeatFlux               = 0.0_RKIND
          surfaceConductiveFlux         = 0.0_RKIND
          surfaceHeatFluxCategory       = 0.0_RKIND
          surfaceConductiveFluxCategory = 0.0_RKIND
          latentHeatFluxCategory        = 0.0_RKIND
          sensibleHeatFluxCategory      = 0.0_RKIND

          ! melt growth rates
          call MPAS_pool_get_subpool(block % structs, "melt_growth_rates", meltGrowthRatesPool)

          call MPAS_pool_get_array(meltGrowthRatesPool, "congelation", congelation)
          call MPAS_pool_get_array(meltGrowthRatesPool, "frazilFormation", frazilFormation)
          call MPAS_pool_get_array(meltGrowthRatesPool, "snowiceFormation", snowiceFormation)
          call MPAS_pool_get_array(meltGrowthRatesPool, "snowThicknessChange", snowThicknessChange)
          call MPAS_pool_get_array(meltGrowthRatesPool, "surfaceIceMelt", surfaceIceMelt)
          call MPAS_pool_get_array(meltGrowthRatesPool, "snowMelt", snowMelt)
          call MPAS_pool_get_array(meltGrowthRatesPool, "basalIceMelt", basalIceMelt)
          call MPAS_pool_get_array(meltGrowthRatesPool, "lateralIceMelt", lateralIceMelt)

          congelation         = 0.0_RKIND
          frazilFormation     = 0.0_RKIND
          snowiceFormation    = 0.0_RKIND
          snowThicknessChange = 0.0_RKIND
          surfaceIceMelt      = 0.0_RKIND
          snowMelt            = 0.0_RKIND
          basalIceMelt        = 0.0_RKIND
          lateralIceMelt      = 0.0_RKIND

          ! tendancies
          call MPAS_pool_get_config(block % configs, "config_use_ice_age", config_use_ice_age)

          call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnosticsPool)
          call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)

          call MPAS_pool_get_array(diagnosticsPool, "iceAreaTendencyThermodynamics", iceAreaTendencyThermodynamics)
          call MPAS_pool_get_array(diagnosticsPool, "iceVolumeTendencyThermodynamics", iceVolumeTendencyThermodynamics)
          call MPAS_pool_get_array(diagnosticsPool, "iceAgeTendencyThermodynamics", iceAgeTendencyThermodynamics)

          call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
          call MPAS_pool_get_array(tracersAggregatePool, "iceVolumeCell", iceVolumeCell)
          call MPAS_pool_get_array(tracersAggregatePool, "iceAgeCell", iceAgeCell)

          ! thermodynamic tendencies
          iceAreaTendencyThermodynamics   = iceAreaCell
          iceVolumeTendencyThermodynamics = iceVolumeCell
          if (config_use_ice_age) then
             iceAgeTendencyThermodynamics = iceAgeCell
          else
             iceAgeTendencyThermodynamics = 0.0_RKIND
          endif

          ! ponds
          call MPAS_pool_get_subpool(block % structs, "ponds", pondsPool)

          call MPAS_pool_get_array(pondsPool, "pondFreshWaterFlux", pondFreshWaterFlux)

          pondFreshWaterFlux(:) = 0.0_RKIND

          !fresh_ai  (:,:,:) = 0.0_RKIND
          !fsalt_ai  (:,:,:) = 0.0_RKIND
          !fhocn_ai  (:,:,:) = 0.0_RKIND
          !fswthru_ai(:,:,:) = 0.0_RKIND

          ! shortwave
          call MPAS_pool_get_subpool(block % structs, "shortwave", shortwavePool)

          call MPAS_pool_get_array(shortwavePool, "bareIceAlbedoCell", bareIceAlbedoCell)
          call MPAS_pool_get_array(shortwavePool, "snowAlbedoCell", snowAlbedoCell)
          call MPAS_pool_get_array(shortwavePool, "pondAlbedoCell", pondAlbedoCell)

          bareIceAlbedoCell = 0.0_RKIND
          snowAlbedoCell    = 0.0_RKIND
          pondAlbedoCell    = 0.0_RKIND

          ! form drag
          call MPAS_pool_get_config(block % configs, "config_ice_ocean_drag_coefficient", &
                                    config_ice_ocean_drag_coefficient)

          call MPAS_pool_get_subpool(block % structs, "drag", dragPool)

          call MPAS_pool_get_array(dragPool, "oceanDragCoefficient", oceanDragCoefficient)
          call MPAS_pool_get_array(dragPool, "airDragCoefficient", airDragCoefficient)

          oceanDragCoefficient = config_ice_ocean_drag_coefficient
          airDragCoefficient   = seaice_icepack_initial_air_drag_coefficient()

          call MPAS_pool_get_config(block % configs, "config_use_form_drag", config_use_form_drag)

          if (config_use_form_drag) then

             call MPAS_pool_get_array(dragPool, "airOceanDragCoefficientRatio", airOceanDragCoefficientRatio)
             call MPAS_pool_get_array(dragPool, "oceanDragCoefficientSkin", oceanDragCoefficientSkin)
             call MPAS_pool_get_array(dragPool, "oceanDragCoefficientFloe", oceanDragCoefficientFloe)
             call MPAS_pool_get_array(dragPool, "oceanDragCoefficientKeel", oceanDragCoefficientKeel)
             call MPAS_pool_get_array(dragPool, "airDragCoefficientSkin", airDragCoefficientSkin)
             call MPAS_pool_get_array(dragPool, "airDragCoefficientFloe", airDragCoefficientFloe)
             call MPAS_pool_get_array(dragPool, "airDragCoefficientPond", airDragCoefficientPond)
             call MPAS_pool_get_array(dragPool, "airDragCoefficientRidge", airDragCoefficientRidge)
             call MPAS_pool_get_array(dragPool, "dragFreeboard", dragFreeboard)
             call MPAS_pool_get_array(dragPool, "dragIceSnowDraft", dragIceSnowDraft)
             call MPAS_pool_get_array(dragPool, "dragRidgeHeight", dragRidgeHeight)
             call MPAS_pool_get_array(dragPool, "dragRidgeSeparation", dragRidgeSeparation)
             call MPAS_pool_get_array(dragPool, "dragKeelDepth", dragKeelDepth)
             call MPAS_pool_get_array(dragPool, "dragKeelSeparation", dragKeelSeparation)
             call MPAS_pool_get_array(dragPool, "dragFloeLength", dragFloeLength)
             call MPAS_pool_get_array(dragPool, "dragFloeSeparation", dragFloeSeparation)

             airOceanDragCoefficientRatio = 0.0_RKIND
             oceanDragCoefficientSkin     = 0.0_RKIND
             oceanDragCoefficientFloe     = 0.0_RKIND
             oceanDragCoefficientKeel     = 0.0_RKIND
             airDragCoefficientSkin       = 0.0_RKIND
             airDragCoefficientFloe       = 0.0_RKIND
             airDragCoefficientPond       = 0.0_RKIND
             airDragCoefficientRidge      = 0.0_RKIND
             dragFreeboard                = 0.0_RKIND
             dragIceSnowDraft             = 0.0_RKIND
             dragRidgeHeight              = 0.0_RKIND
             dragRidgeSeparation          = 0.0_RKIND
             dragKeelDepth                = 0.0_RKIND
             dragKeelSeparation           = 0.0_RKIND
             dragFloeLength               = 0.0_RKIND
             dragFloeSeparation           = 0.0_RKIND

          endif ! config_use_form_drag

          ! snow
          call MPAS_pool_get_config(block % configs, "config_use_column_snow_tracers", config_use_column_snow_tracers)
          if (config_use_column_snow_tracers) then

             call MPAS_pool_get_subpool(block % structs, "snow", snowPool)
             call MPAS_pool_get_array(snowPool, "snowLossToLeads", snowLossToLeads)
             call MPAS_pool_get_array(snowPool, "snowMeltMassCell", snowMeltMassCell)
             call MPAS_pool_get_array(snowPool, "snowMeltMassCategory", snowMeltMassCategory)
             call MPAS_pool_get_array(snowPool, "snowDensityViaContent", snowDensityViaContent)
             call MPAS_pool_get_array(snowPool, "snowDensityViaCompaction", snowDensityViaCompaction)
             call MPAS_pool_get_array(snowPool, "snowRadiusInStandardRadiationScheme", snowRadiusInStandardRadiationScheme)
             call MPAS_pool_get_array(snowPool, "snowRadiusInStandardRadiationSchemeCategory", snowRadiusInStandardRadiationSchemeCategory)

             snowLossToLeads             = 0.0_RKIND
             snowMeltMassCell            = 0.0_RKIND
             snowMeltMassCategory        = 0.0_RKIND
             snowDensityViaContent       = 0.0_RKIND
             snowDensityViaCompaction    = 0.0_RKIND
             snowRadiusInStandardRadiationScheme         = 0.0_RKIND
             snowRadiusInStandardRadiationSchemeCategory = 0.0_RKIND

          end if ! config_use_column_snow_tracers

          block => block % next
       end do

    endif ! config_use_column_physics

  end subroutine seaice_icepack_reinitialize_diagnostics_thermodynamics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack_reinitialize_diagnostics_dynamics
!
!> \brief Reinitialize dynamics diagnostics
!> \author Adrian K. Turner, LANL
!> \date 27th September 2015
!> \details
!>  Reinitialize dynamics diagnostics
!
!-----------------------------------------------------------------------

  subroutine seaice_icepack_reinitialize_diagnostics_dynamics(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         velocityWeakPool, &
         velocityVariationalPool, &
         ridgingPool, &
         diagnosticsPool, &
         tracersAggregatePool

    ! dynamics
    real(kind=RKIND), dimension(:), pointer :: &
         oceanStressU, &
         oceanStressV, &
         airStressVertexU, &
         airStressVertexV, &
         stressDivergenceU, &
         stressDivergenceV, &
         surfaceTiltForceU, &
         surfaceTiltForceV

    real(kind=RKIND), dimension(:,:), pointer :: &
         principalStress1Var, &
         principalStress2Var, &
         replacementPressureVar

    real(kind=RKIND), dimension(:), pointer :: &
         principalStress1Weak, &
         principalStress2Weak, &
         replacementPressureWeak

    ! ridging
    real(kind=RKIND), dimension(:), pointer :: &
         areaLossRidge, &
         areaGainRidge, &
         iceVolumeRidged, &
         openingRateRidge

    ! diagnostic tendencies
    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaTendencyTransport, &
         iceVolumeTendencyTransport, &
         iceAgeTendencyTransport, &
         iceAreaCell, &
         iceVolumeCell, &
         iceAgeCell

    logical, pointer :: &
         config_use_ice_age, &
         config_use_column_physics, &
         config_use_velocity_solver

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_physics", config_use_column_physics)
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_velocity_solver", config_use_velocity_solver)
    call MPAS_pool_get_config(domain % blocklist % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)

    if (config_use_column_physics) then

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

          call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
          call MPAS_pool_get_array(velocitySolverPool, "airStressVertexU", airStressVertexU)
          call MPAS_pool_get_array(velocitySolverPool, "airStressVertexV", airStressVertexV)
          call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
          call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)
          call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
          call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)

          oceanStressU        = 0.0_RKIND
          oceanStressV        = 0.0_RKIND
          airStressVertexU    = 0.0_RKIND
          airStressVertexV    = 0.0_RKIND
          stressDivergenceU   = 0.0_RKIND
          stressDivergenceV   = 0.0_RKIND
          surfaceTiltForceU   = 0.0_RKIND
          surfaceTiltForceV   = 0.0_RKIND

          if (config_use_velocity_solver .and. trim(config_stress_divergence_scheme) == "weak") then

             call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocityWeakPool)

             call MPAS_pool_get_array(velocityWeakPool, "principalStress1", principalStress1Weak)
             call MPAS_pool_get_array(velocityWeakPool, "principalStress2", principalStress2Weak)
             call MPAS_pool_get_array(velocityWeakPool, "replacementPressure", replacementPressureWeak)

             principalStress1Weak    = 0.0_RKIND
             principalStress2Weak    = 0.0_RKIND
             replacementPressureWeak = 0.0_RKIND

          else if (config_use_velocity_solver .and. trim(config_stress_divergence_scheme) == "variational") then

             call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocityVariationalPool)

             call MPAS_pool_get_array(velocityVariationalPool, "principalStress1", principalStress1Var)
             call MPAS_pool_get_array(velocityVariationalPool, "principalStress2", principalStress2Var)
             call MPAS_pool_get_array(velocityVariationalPool, "replacementPressure", replacementPressureVar)

             principalStress1Var    = 0.0_RKIND
             principalStress2Var    = 0.0_RKIND
             replacementPressureVar = 0.0_RKIND

          endif

          call MPAS_pool_get_subpool(block % structs, "ridging", ridgingPool)

          call MPAS_pool_get_array(ridgingPool, "areaLossRidge", areaLossRidge)
          call MPAS_pool_get_array(ridgingPool, "areaGainRidge", areaGainRidge)
          call MPAS_pool_get_array(ridgingPool, "iceVolumeRidged", iceVolumeRidged)
          call MPAS_pool_get_array(ridgingPool, "openingRateRidge", openingRateRidge)

          areaLossRidge    = 0.0_RKIND
          areaGainRidge    = 0.0_RKIND
          iceVolumeRidged  = 0.0_RKIND
          openingRateRidge = 0.0_RKIND

          ! tendancies
          call MPAS_pool_get_config(block % configs, "config_use_ice_age", config_use_ice_age)

          call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnosticsPool)
          call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)

          call MPAS_pool_get_array(diagnosticsPool, "iceAreaTendencyTransport", iceAreaTendencyTransport)
          call MPAS_pool_get_array(diagnosticsPool, "iceVolumeTendencyTransport", iceVolumeTendencyTransport)
          call MPAS_pool_get_array(diagnosticsPool, "iceAgeTendencyTransport", iceAgeTendencyTransport)

          call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
          call MPAS_pool_get_array(tracersAggregatePool, "iceVolumeCell", iceVolumeCell)
          call MPAS_pool_get_array(tracersAggregatePool, "iceAgeCell", iceAgeCell)

          ! transport tendencies
          iceAreaTendencyTransport   = iceAreaCell
          iceVolumeTendencyTransport = iceVolumeCell
          if (config_use_ice_age) then
             iceAgeTendencyTransport = iceAgeCell
          else
             iceAgeTendencyTransport = 0.0_RKIND
          endif

          block => block % next
       end do

    endif ! config_use_column_physics

  end subroutine seaice_icepack_reinitialize_diagnostics_dynamics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_icepack_reinitialize_diagnostics_bgc
!
!> \brief Reinitialize BGC diagnostics
!> \author Adrian K. Turner, Nicole Jeffery, LANL
!> \date 2022-2023
!> \details
!>  Reinitialize BGC diagnostics
!
!-----------------------------------------------------------------------

  subroutine seaice_icepack_reinitialize_diagnostics_bgc(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         biogeochemistryPool, &
         diagnostics_biogeochemistryPool, &
         tracers_aggregatePool

    ! biogeochemistry
    real(kind=RKIND), dimension(:), pointer :: &
         primaryProduction, &
         netSpecificAlgalGrowthRate, &
         netBrineHeight, &
!         zSalinityFlux, &        !echmod deprecate
!         zSalinityGDFlux, &      !echmod deprecate
         totalChlorophyll, &
         totalCarbonContentCell, &
         totalVerticalDiatomIce, &
         totalVerticalSmallPlanktonIce, &
         totalVerticalPhaeocystisIce, &
         totalVerticalNitrateIce, &
         totalVerticalPolysaccharidsIce, &
         totalVerticalLipidsIce, &
         totalVerticalDICIce, &
         totalVerticalAmmoniumIce, &
         totalVerticalSilicateIce, &
         totalVerticalDMSPpIce, &
         totalVerticalDMSPdIce, &
         totalVerticalDMSIce, &
         totalVerticalNonreactiveIce, &
         totalVerticalHumicsIce, &
         totalVerticalProteinsIce, &
         totalVerticalDissolvedIronIce, &
         totalVerticalParticulateIronIce, &
         totalVerticalDustIce, &
         totalVerticalDustSnow, &
         totalVerticalBC1Ice, &
         totalVerticalBC1Snow, &
         totalVerticalBC2Ice, &
         totalVerticalBC2Snow, &
         totalVerticalBCIce, &
         totalVerticalBCSnow, &
         totalVerticalAlgaeCarbonIce, &
         totalVerticalDOCLabileIce, &
         totalVerticalDissolvedIronSnow

    real(kind=RKIND), dimension(:,:), pointer :: &
         oceanBioFluxes, &
         atmosIceBioFluxes, &
         snowIceBioFluxes, &
         totalVerticalBiologyIce, &
         totalVerticalBiologySnow, &
         bgridPorosityIceCell, &
         bgridSalinityIceCell, &
         bgridTemperatureIceCell, &
         verticalBCTotalIceCell, &
         verticalDustTotalIceCell, &
         verticalBCTotalSnowCell, &
         verticalDustTotalSnowCell


    real(kind=RKIND), dimension(:,:,:), pointer :: &
         bioTracerShortwave

    logical, pointer :: &
         config_use_column_biogeochemistry, &
         config_use_column_shortwave, &
         config_use_column_physics, &
         config_use_vertical_tracers, &
!         config_use_vertical_zsalinity, &
         config_use_zaerosols

    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_physics", config_use_column_physics)

     if (config_use_column_physics) then

       block => domain % blocklist
       do while (associated(block))

          ! biogeochemistry
          call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
          call MPAS_pool_get_config(block % configs, "config_use_zaerosols", config_use_zaerosols)
          call MPAS_pool_get_config(block % configs, "config_use_vertical_tracers", config_use_vertical_tracers)
!          call MPAS_pool_get_config(block % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)

          if (config_use_column_biogeochemistry .or. config_use_zaerosols) then

             call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistryPool)
             call MPAS_pool_get_subpool(block % structs, "diagnostics_biogeochemistry", diagnostics_biogeochemistryPool)
             call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregatePool)

             if (config_use_vertical_tracers) then
                call MPAS_pool_get_array(biogeochemistryPool, "primaryProduction", primaryProduction)
                call MPAS_pool_get_array(biogeochemistryPool, "totalChlorophyll", totalChlorophyll)
                call MPAS_pool_get_array(biogeochemistryPool, "netSpecificAlgalGrowthRate", netSpecificAlgalGrowthRate)
                call MPAS_pool_get_array(diagnostics_biogeochemistryPool, "bgridSalinityIceCell", bgridSalinityIceCell)
                call MPAS_pool_get_array(diagnostics_biogeochemistryPool, "bgridPorosityIceCell", bgridPorosityIceCell)
                call MPAS_pool_get_array(diagnostics_biogeochemistryPool, "bgridTemperatureIceCell", bgridTemperatureIceCell)
                primaryProduction             = 0.0_RKIND
                totalChlorophyll              = 0.0_RKIND
                netSpecificAlgalGrowthRate    = 0.0_RKIND
                bgridSalinityIceCell          = 0.0_RKIND
                bgridPorosityIceCell          = 0.0_RKIND
                bgridTemperatureIceCell       = 0.0_RKIND

             end if

             call MPAS_pool_get_array(biogeochemistryPool, "netBrineHeight", netBrineHeight)
             call MPAS_pool_get_array(biogeochemistryPool, "oceanBioFluxes", oceanBioFluxes)
             call MPAS_pool_get_array(biogeochemistryPool, "atmosIceBioFluxes", atmosIceBioFluxes)
             call MPAS_pool_get_array(biogeochemistryPool, "snowIceBioFluxes", snowIceBioFluxes)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalBiologyIce", totalVerticalBiologyIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalBiologySnow", totalVerticalBiologySnow)
             call MPAS_pool_get_array(biogeochemistryPool, "totalCarbonContentCell", totalCarbonContentCell)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalDiatomIce", totalVerticalDiatomIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalSmallPlanktonIce", totalVerticalSmallPlanktonIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalPhaeocystisIce", totalVerticalPhaeocystisIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalPolysaccharidsIce", totalVerticalPolysaccharidsIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalLipidsIce", totalVerticalLipidsIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalDICIce", totalVerticalDICIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalProteinsIce", totalVerticalProteinsIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalNitrateIce", totalVerticalNitrateIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalSilicateIce", totalVerticalSilicateIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalAmmoniumIce", totalVerticalAmmoniumIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalDMSIce", totalVerticalDMSIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalDMSPpIce", totalVerticalDMSPpIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalDMSPdIce", totalVerticalDMSPdIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalNonreactiveIce", totalVerticalNonreactiveIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalHumicsIce", totalVerticalHumicsIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalParticulateIronIce", totalVerticalParticulateIronIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalDissolvedIronIce", totalVerticalDissolvedIronIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalDissolvedIronSnow", totalVerticalDissolvedIronSnow)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalDustIce", totalVerticalDustIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalDustSnow", totalVerticalDustSnow)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalBC1Ice", totalVerticalBC1Ice)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalBC1Snow", totalVerticalBC1Snow)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalBC2Ice", totalVerticalBC2Ice)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalBC2Snow", totalVerticalBC2Snow)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalBCIce", totalVerticalBCIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalBCSnow", totalVerticalBCSnow)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalAlgaeCarbonIce", totalVerticalAlgaeCarbonIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalDOCLabileIce", totalVerticalDOCLabileIce)
             call MPAS_pool_get_array(tracers_aggregatePool, "verticalBCTotalIceCell", verticalBCTotalIceCell)
             call MPAS_pool_get_array(tracers_aggregatePool, "verticalDustTotalIceCell", verticalDustTotalIceCell)
             call MPAS_pool_get_array(tracers_aggregatePool, "verticalDustTotalSnowCell", verticalDustTotalSnowCell)
             call MPAS_pool_get_array(tracers_aggregatePool, "verticalBCTotalSnowCell", verticalBCTotalSnowCell)


             netBrineHeight             = 0.0_RKIND
             oceanBioFluxes             = 0.0_RKIND
             atmosIceBioFluxes          = 0.0_RKIND
             snowIceBioFluxes           = 0.0_RKIND
             totalVerticalBiologyIce    = 0.0_RKIND
             totalVerticalBiologySnow   = 0.0_RKIND
             totalCarbonContentCell     = 0.0_RKIND
             totalVerticalDiatomIce       = 0.0_RKIND
             totalVerticalSmallPlanktonIce = 0.0_RKIND
             totalVerticalPhaeocystisIce   = 0.0_RKIND
             totalVerticalPolysaccharidsIce = 0.0_RKIND
             totalVerticalLipidsIce        = 0.0_RKIND
             totalVerticalDICIce           = 0.0_RKIND
             totalVerticalProteinsIce      = 0.0_RKIND
             totalVerticalNitrateIce       = 0.0_RKIND
             totalVerticalSilicateIce      = 0.0_RKIND
             totalVerticalAmmoniumIce      = 0.0_RKIND
             totalVerticalDMSIce           = 0.0_RKIND
             totalVerticalDMSPpIce         = 0.0_RKIND
             totalVerticalDMSPdIce         = 0.0_RKIND
             totalVerticalNonreactiveIce   = 0.0_RKIND
             totalVerticalHumicsIce        = 0.0_RKIND
             totalVerticalParticulateIronIce = 0.0_RKIND
             totalVerticalDissolvedIronIce    = 0.0_RKIND
             totalVerticalDissolvedIronSnow   = 0.0_RKIND
             totalVerticalBC1Ice              = 0.0_RKIND
             totalVerticalBC2Ice           = 0.0_RKIND
             totalVerticalDustIce          = 0.0_RKIND
             totalVerticalBC1Snow          = 0.0_RKIND
             totalVerticalBC2Snow          = 0.0_RKIND
             totalVerticalBCSnow           = 0.0_RKIND
             totalVerticalBCIce            = 0.0_RKIND
             totalVerticalDustSnow         = 0.0_RKIND
             totalVerticalAlgaeCarbonIce   = 0.0_RKIND
             totalVerticalDOCLabileIce     = 0.0_RKIND
             verticalBCTotalIceCell        = 0.0_RKIND
             verticalDustTotalIceCell      = 0.0_RKIND
             verticalBCTotalSnowCell       = 0.0_RKIND
             verticalDustTotalSnowCell     = 0.0_RKIND

          endif

          call MPAS_pool_get_config(block % configs, "config_use_column_shortwave", config_use_column_shortwave)

          if (config_use_column_biogeochemistry .or. config_use_column_shortwave .or. config_use_zaerosols) then

             call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistryPool)
             call MPAS_pool_get_array(biogeochemistryPool, "bioTracerShortwave", bioTracerShortwave)
             bioTracerShortwave         = 0.0_RKIND

          endif

          block => block % next
       end do

    endif ! config_use_column_physics

  end subroutine seaice_icepack_reinitialize_diagnostics_bgc

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_separate_snow_ice_tracers
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 13 March 2017
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine  column_separate_snow_ice_tracers(domain)

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers

    logical, pointer :: &
         config_use_vertical_biochemistry, &
         config_use_nitrate, &
         config_use_carbon, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_nonreactive, &
         config_use_humics, &
         config_use_DON, &
         config_use_iron, &
         config_use_zaerosols

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         verticalAerosolsConc, &
         verticalAerosolsSnow, &
         verticalAerosolsIce, &
         verticalDissolvedIronConc, &
         verticalParticulateIronConc, &
         verticalHumicsConc, &
         verticalNonreactiveConc, &
         verticalDMSPdConc, &
         verticalDMSPpConc, &
         verticalDMSConc, &
         verticalAmmoniumConc, &
         verticalSilicateConc, &
         verticalNitrateConc, &
         verticalDONConc, &
         verticalDICConc, &
         verticalDOCConc, &
         verticalAlgaeConc, &
         verticalDissolvedIronSnow, &
         verticalParticulateIronSnow, &
         verticalHumicsSnow, &
         verticalNonreactiveSnow, &
         verticalDMSPdSnow, &
         verticalDMSPpSnow, &
         verticalDMSSnow, &
         verticalAmmoniumSnow, &
         verticalSilicateSnow, &
         verticalNitrateSnow, &
         verticalDONSnow, &
         verticalDICSnow, &
         verticalDOCSnow, &
         verticalAlgaeSnow, &
         verticalDissolvedIronIce, &
         verticalParticulateIronIce, &
         verticalHumicsIce, &
         verticalNonreactiveIce, &
         verticalDMSPdIce, &
         verticalDMSPpIce, &
         verticalDMSIce, &
         verticalAmmoniumIce, &
         verticalSilicateIce, &
         verticalNitrateIce, &
         verticalDONIce, &
         verticalDICIce, &
         verticalDOCIce, &
         verticalAlgaeIce

    integer, pointer :: &
         nCellsSolve, &
         nBioLayersP1, &
         nBioLayersP3, &
         nCategories, &
         nzAerosols, &
         TWO, &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON, &
         nParticulateIron, &
         nDissolvedIron

    integer :: &
         iCell, &
         iBioTracers, &
         iCategory, &
         iBioData, &
         iBioCount, &
         iSnowCount, &
         iIceCount

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nCategories", nCategories)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
       call MPAS_pool_get_config(block % configs, "config_use_nitrate", config_use_nitrate)
       call MPAS_pool_get_config(block % configs, "config_use_carbon", config_use_carbon)
       call MPAS_pool_get_config(block % configs, "config_use_ammonium",config_use_ammonium)
       call MPAS_pool_get_config(block % configs, "config_use_silicate",config_use_silicate)
       call MPAS_pool_get_config(block % configs, "config_use_DMS",config_use_DMS)
       call MPAS_pool_get_config(block % configs, "config_use_nonreactive",config_use_nonreactive)
       call MPAS_pool_get_config(block % configs, "config_use_humics",config_use_humics)
       call MPAS_pool_get_config(block % configs, "config_use_DON",config_use_DON)
       call MPAS_pool_get_config(block % configs, "config_use_iron",config_use_iron)
       call MPAS_pool_get_config(block % configs, "config_use_zaerosols",config_use_zaerosols)

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nzAerosols", nzAerosols)
       call MPAS_pool_get_dimension(mesh, "nBioLayersP1", nBioLayersP1)
       call MPAS_pool_get_dimension(mesh, "nBioLayersP3", nBioLayersP3)
       call MPAS_pool_get_dimension(mesh, "TWO", TWO)
       call MPAS_pool_get_dimension(mesh, "nAlgae", nAlgae)
       call MPAS_pool_get_dimension(mesh, "nDOC", nDOC)
       call MPAS_pool_get_dimension(mesh, "nDIC", nDIC)
       call MPAS_pool_get_dimension(mesh, "nDON", nDON)
       call MPAS_pool_get_dimension(mesh, "nParticulateIron", nParticulateIron)
       call MPAS_pool_get_dimension(mesh, "nDissolvedIron", nDissolvedIron)

       call MPAS_pool_get_array(tracers, "verticalAerosolsConc",verticalAerosolsConc,1)
       call MPAS_pool_get_array(tracers, "verticalAerosolsSnow",verticalAerosolsSnow,1)
       call MPAS_pool_get_array(tracers, "verticalAerosolsIce",verticalAerosolsIce,1)
       call MPAS_pool_get_array(tracers, "verticalAlgaeConc",verticalAlgaeConc,1)
       call MPAS_pool_get_array(tracers, "verticalDOCConc",verticalDOCConc,1)
       call MPAS_pool_get_array(tracers, "verticalDICConc",verticalDICConc,1)
       call MPAS_pool_get_array(tracers, "verticalDONConc",verticalDONConc,1)
       call MPAS_pool_get_array(tracers, "verticalNitrateConc",verticalNitrateConc,1)
       call MPAS_pool_get_array(tracers, "verticalSilicateConc",verticalSilicateConc,1)
       call MPAS_pool_get_array(tracers, "verticalAmmoniumConc",verticalAmmoniumConc,1)
       call MPAS_pool_get_array(tracers, "verticalDMSConc",verticalDMSConc,1)
       call MPAS_pool_get_array(tracers, "verticalDMSPpConc",verticalDMSPpConc,1)
       call MPAS_pool_get_array(tracers, "verticalDMSPdConc",verticalDMSPdConc,1)
       call MPAS_pool_get_array(tracers, "verticalNonreactiveConc",verticalNonreactiveConc,1)
       call MPAS_pool_get_array(tracers, "verticalHumicsConc",verticalHumicsConc,1)
       call MPAS_pool_get_array(tracers, "verticalParticulateIronConc",verticalParticulateIronConc,1)
       call MPAS_pool_get_array(tracers, "verticalDissolvedIronConc",verticalDissolvedIronConc,1)
       call MPAS_pool_get_array(tracers, "verticalAlgaeSnow",verticalAlgaeSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDOCSnow",verticalDOCSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDICSnow",verticalDICSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDONSnow",verticalDONSnow,1)
       call MPAS_pool_get_array(tracers, "verticalNitrateSnow",verticalNitrateSnow,1)
       call MPAS_pool_get_array(tracers, "verticalSilicateSnow",verticalSilicateSnow,1)
       call MPAS_pool_get_array(tracers, "verticalAmmoniumSnow",verticalAmmoniumSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDMSSnow",verticalDMSSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDMSPpSnow",verticalDMSPpSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDMSPdSnow",verticalDMSPdSnow,1)
       call MPAS_pool_get_array(tracers, "verticalNonreactiveSnow",verticalNonreactiveSnow,1)
       call MPAS_pool_get_array(tracers, "verticalHumicsSnow",verticalHumicsSnow,1)
       call MPAS_pool_get_array(tracers, "verticalParticulateIronSnow",verticalParticulateIronSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDissolvedIronSnow",verticalDissolvedIronSnow,1)
       call MPAS_pool_get_array(tracers, "verticalAlgaeIce",verticalAlgaeIce,1)
       call MPAS_pool_get_array(tracers, "verticalDOCIce",verticalDOCIce,1)
       call MPAS_pool_get_array(tracers, "verticalDICIce",verticalDICIce,1)
       call MPAS_pool_get_array(tracers, "verticalDONIce",verticalDONIce,1)
       call MPAS_pool_get_array(tracers, "verticalNitrateIce",verticalNitrateIce,1)
       call MPAS_pool_get_array(tracers, "verticalSilicateIce",verticalSilicateIce,1)
       call MPAS_pool_get_array(tracers, "verticalAmmoniumIce",verticalAmmoniumIce,1)
       call MPAS_pool_get_array(tracers, "verticalDMSIce",verticalDMSIce,1)
       call MPAS_pool_get_array(tracers, "verticalDMSPpIce",verticalDMSPpIce,1)
       call MPAS_pool_get_array(tracers, "verticalDMSPdIce",verticalDMSPdIce,1)
       call MPAS_pool_get_array(tracers, "verticalNonreactiveIce",verticalNonreactiveIce,1)
       call MPAS_pool_get_array(tracers, "verticalHumicsIce",verticalHumicsIce,1)
       call MPAS_pool_get_array(tracers, "verticalParticulateIronIce",verticalParticulateIronIce,1)
       call MPAS_pool_get_array(tracers, "verticalDissolvedIronIce",verticalDissolvedIronIce,1)

       do iCell = 1, nCellsSolve
          do iCategory = 1, nCategories

             ! aerosols
             if (config_use_zaerosols) then
                do iBioTracers = 1, nzAerosols
                   iSnowCount = (iBioTracers-1)*2
                   iIceCount = (iBioTracers-1)*nBioLayersP1
                   iBioData =  (iBioTracers-1)*nBioLayersP3

                   do iBioCount = 1,TWO
                      verticalAerosolsSnow(iBioCount+iSnowCount,iCategory,iCell) = &
                         verticalAerosolsConc(iBioData+nBioLayersP1+iBioCount,iCategory,iCell)

                   enddo
                   do iBioCount = 1, nBioLayersP1
                      verticalAerosolsIce(iBioCount+iIceCount,iCategory,iCell) = &
                         verticalAerosolsConc(iBioData+iBioCount,iCategory,iCell)
                   enddo
                enddo
             endif

             ! algal nitrogen
             if (config_use_vertical_biochemistry) then
                do iBioTracers = 1, nAlgae
                   iSnowCount = (iBioTracers-1)*2
                   iIceCount = (iBioTracers-1)*nBioLayersP1
                   iBioData =  (iBioTracers-1)*nBioLayersP3

                   do iBioCount = 1,TWO
                      verticalAlgaeSnow(iBioCount+iSnowCount,iCategory,iCell) = &
                         verticalAlgaeConc(iBioData+nBioLayersP1+iBioCount,iCategory,iCell)
                   enddo
                   do iBioCount = 1, nBioLayersP1
                      verticalAlgaeIce(iBioCount+iIceCount,iCategory,iCell) = &
                         verticalAlgaeConc(iBioData+iBioCount,iCategory,iCell)
                   enddo
                enddo
             endif

             ! dissolved organic and inorganic carbon
             if (config_use_carbon) then
                do iBioTracers = 1, nDOC
                   iSnowCount = (iBioTracers-1)*2
                   iIceCount = (iBioTracers-1)*nBioLayersP1
                   iBioData =  (iBioTracers-1)*nBioLayersP3

                   do iBioCount = 1,TWO
                      verticalDOCSnow(iBioCount+iSnowCount,iCategory,iCell) = &
                         verticalDOCConc(iBioData+nBioLayersP1+iBioCount,iCategory,iCell)
                   enddo
                   do iBioCount = 1, nBioLayersP1
                      verticalDOCIce(iBioCount+iIceCount,iCategory,iCell) = &
                         verticalDOCConc(iBioData+iBioCount,iCategory,iCell)
                   enddo
                enddo
                do iBioTracers = 1, nDIC
                   iSnowCount = (iBioTracers-1)*2
                   iIceCount = (iBioTracers-1)*nBioLayersP1
                   iBioData =  (iBioTracers-1)*nBioLayersP3

                   do iBioCount = 1,TWO
                      verticalDICSnow(iBioCount+iSnowCount,iCategory,iCell) = &
                         verticalDICConc(iBioData+nBioLayersP1+iBioCount,iCategory,iCell)
                   enddo
                   do iBioCount = 1, nBioLayersP1
                      verticalDICIce(iBioCount+iIceCount,iCategory,iCell) = &
                         verticalDICConc(iBioData+iBioCount,iCategory,iCell)
                   enddo
                enddo
             endif

             ! nitrate
             if (config_use_nitrate) then
                do iBioCount = 1,TWO
                   verticalNitrateSnow(iBioCount,iCategory,iCell) = &
                      verticalNitrateConc(nBioLayersP1+iBioCount,iCategory,iCell)
                enddo
                do iBioCount = 1, nBioLayersP1
                   verticalNitrateIce(iBioCount,iCategory,iCell) = &
                      verticalNitrateConc(iBioCount,iCategory,iCell)
                enddo
             endif

             ! ammonium
             if (config_use_ammonium) then
                do iBioCount = 1,TWO
                   verticalAmmoniumSnow(iBioCount,iCategory,iCell) = &
                      verticalAmmoniumConc(nBioLayersP1+iBioCount,iCategory,iCell)
                enddo
                do iBioCount = 1, nBioLayersP1
                   verticalAmmoniumIce(iBioCount,iCategory,iCell) = &
                      verticalAmmoniumConc(iBioCount,iCategory,iCell)
                enddo
             endif

             ! silicate
             if (config_use_silicate) then
                do iBioCount = 1,TWO
                   verticalSilicateSnow(iBioCount,iCategory,iCell) = &
                      verticalSilicateConc(nBioLayersP1+iBioCount,iCategory,iCell)
                enddo
                do iBioCount = 1, nBioLayersP1
                   verticalSilicateIce(iBioCount,iCategory,iCell) = &
                      verticalSilicateConc(iBioCount,iCategory,iCell)
                enddo
             endif

             ! DMS, DMSPp, DMPSd
             if (config_use_DMS) then
                do iBioCount = 1,TWO
                   verticalDMSSnow(iBioCount,iCategory,iCell) = &
                      verticalDMSConc(nBioLayersP1+iBioCount,iCategory,iCell)
                   verticalDMSPpSnow(iBioCount,iCategory,iCell) = &
                      verticalDMSPpConc(nBioLayersP1+iBioCount,iCategory,iCell)
                   verticalDMSPdSnow(iBioCount,iCategory,iCell) = &
                      verticalDMSPdConc(nBioLayersP1+iBioCount,iCategory,iCell)
                enddo
                do iBioCount = 1, nBioLayersP1
                   verticalDMSIce(iBioCount,iCategory,iCell) = &
                      verticalDMSConc(iBioCount,iCategory,iCell)
                   verticalDMSPpIce(iBioCount,iCategory,iCell) = &
                      verticalDMSPpConc(iBioCount,iCategory,iCell)
                   verticalDMSPdIce(iBioCount,iCategory,iCell) = &
                      verticalDMSPdConc(iBioCount,iCategory,iCell)
                enddo
             endif

             ! nonreactive tracer
             if (config_use_nonreactive) then
                do iBioCount = 1,TWO
                   verticalNonreactiveSnow(iBioCount,iCategory,iCell) = &
                      verticalNonreactiveConc(nBioLayersP1+iBioCount,iCategory,iCell)
                enddo
                do iBioCount = 1, nBioLayersP1
                   verticalNonreactiveIce(iBioCount,iCategory,iCell) = &
                      verticalNonreactiveConc(iBioCount,iCategory,iCell)
                enddo
             endif

             ! humics
             if (config_use_humics) then
                do iBioCount = 1,TWO
                   verticalHumicsSnow(iBioCount,iCategory,iCell) = &
                      verticalHumicsConc(nBioLayersP1+iBioCount,iCategory,iCell)
                enddo
                do iBioCount = 1, nBioLayersP1
                   verticalHumicsIce(iBioCount,iCategory,iCell) = &
                      verticalHumicsConc(iBioCount,iCategory,iCell)
                enddo
             endif

             ! proteins and amino acids
             if (config_use_DON) then
                do iBioTracers = 1, nDON
                   iSnowCount = (iBioTracers-1)*2
                   iIceCount = (iBioTracers-1)*nBioLayersP1
                   iBioData =  (iBioTracers-1)*nBioLayersP3

                   do iBioCount = 1,TWO
                      verticalDONSnow(iBioCount+iSnowCount,iCategory,iCell) = &
                         verticalDONConc(iBioData+nBioLayersP1+iBioCount,iCategory,iCell)
                   enddo
                   do iBioCount = 1, nBioLayersP1
                      verticalDONIce(iBioCount+iIceCount,iCategory,iCell) = &
                         verticalDONConc(iBioData+iBioCount,iCategory,iCell)
                   enddo
                enddo
             endif

             ! particulate and dissolved iron
             if (config_use_iron) then
                do iBioTracers = 1, nParticulateIron
                   iSnowCount = (iBioTracers-1)*2
                   iIceCount = (iBioTracers-1)*nBioLayersP1
                   iBioData =  (iBioTracers-1)*nBioLayersP3

                   do iBioCount = 1,TWO
                      verticalParticulateIronSnow(iBioCount+iSnowCount,iCategory,iCell) = &
                         verticalParticulateIronConc(iBioData+nBioLayersP1+iBioCount,iCategory,iCell)
                   enddo
                   do iBioCount = 1, nBioLayersP1
                      verticalParticulateIronIce(iBioCount+iIceCount,iCategory,iCell) = &
                         verticalParticulateIronConc(iBioData+iBioCount,iCategory,iCell)
                   enddo
                enddo
                do iBioTracers = 1, nDissolvedIron
                   iSnowCount = (iBioTracers-1)*2
                   iIceCount = (iBioTracers-1)*nBioLayersP1
                   iBioData =  (iBioTracers-1)*nBioLayersP3

                   do iBioCount = 1,TWO
                      verticalDissolvedIronSnow(iBioCount+iSnowCount,iCategory,iCell) = &
                         verticalDissolvedIronConc(iBioData+nBioLayersP1+iBioCount,iCategory,iCell)
                   enddo
                   do iBioCount = 1, nBioLayersP1
                      verticalDissolvedIronIce(iBioCount+iIceCount,iCategory,iCell) = &
                         verticalDissolvedIronConc(iBioData+iBioCount,iCategory,iCell)
                   enddo
                enddo
             endif

          enddo
       enddo

       block => block % next
    end do

  end subroutine column_separate_snow_ice_tracers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_combine_snow_ice_tracers
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 13 Mar 2017
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine  column_combine_snow_ice_tracers(domain)

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers

    logical, pointer :: &
         config_use_vertical_biochemistry, &
         config_use_nitrate, &
         config_use_carbon, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_nonreactive, &
         config_use_humics, &
         config_use_DON, &
         config_use_iron, &
         config_use_zaerosols

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         verticalAerosolsConc, &
         verticalAerosolsSnow, &
         verticalAerosolsIce, &
         verticalDissolvedIronConc, &
         verticalParticulateIronConc, &
         verticalHumicsConc, &
         verticalNonreactiveConc, &
         verticalDMSPdConc, &
         verticalDMSPpConc, &
         verticalDMSConc, &
         verticalAmmoniumConc, &
         verticalSilicateConc, &
         verticalNitrateConc, &
         verticalDONConc, &
         verticalDICConc, &
         verticalDOCConc, &
         verticalAlgaeConc, &
         verticalDissolvedIronSnow, &
         verticalParticulateIronSnow, &
         verticalHumicsSnow, &
         verticalNonreactiveSnow, &
         verticalDMSPdSnow, &
         verticalDMSPpSnow, &
         verticalDMSSnow, &
         verticalAmmoniumSnow, &
         verticalSilicateSnow, &
         verticalNitrateSnow, &
         verticalDONSnow, &
         verticalDICSnow, &
         verticalDOCSnow, &
         verticalAlgaeSnow, &
         verticalDissolvedIronIce, &
         verticalParticulateIronIce, &
         verticalHumicsIce, &
         verticalNonreactiveIce, &
         verticalDMSPdIce, &
         verticalDMSPpIce, &
         verticalDMSIce, &
         verticalAmmoniumIce, &
         verticalSilicateIce, &
         verticalNitrateIce, &
         verticalDONIce, &
         verticalDICIce, &
         verticalDOCIce, &
         verticalAlgaeIce

    integer, pointer :: &
         nCellsSolve, &
         nBioLayersP1, &
         nBioLayersP3, &
         nCategories, &
         nzAerosols, &
         TWO, &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON, &
         nParticulateIron, &
         nDissolvedIron

    integer :: &
         iCell, &
         iBioTracers, &
         iCategory, &
         iBioData, &
         iBioCount, &
         iSnowCount, &
         iIceCount

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nCategories", nCategories)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
       call MPAS_pool_get_config(block % configs, "config_use_nitrate", config_use_nitrate)
       call MPAS_pool_get_config(block % configs, "config_use_carbon", config_use_carbon)
       call MPAS_pool_get_config(block % configs, "config_use_ammonium",config_use_ammonium)
       call MPAS_pool_get_config(block % configs, "config_use_silicate",config_use_silicate)
       call MPAS_pool_get_config(block % configs, "config_use_DMS",config_use_DMS)
       call MPAS_pool_get_config(block % configs, "config_use_nonreactive",config_use_nonreactive)
       call MPAS_pool_get_config(block % configs, "config_use_humics",config_use_humics)
       call MPAS_pool_get_config(block % configs, "config_use_DON",config_use_DON)
       call MPAS_pool_get_config(block % configs, "config_use_iron",config_use_iron)
       call MPAS_pool_get_config(block % configs, "config_use_zaerosols",config_use_zaerosols)

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nzAerosols", nzAerosols)
       call MPAS_pool_get_dimension(mesh, "nBioLayersP1", nBioLayersP1)
       call MPAS_pool_get_dimension(mesh, "nBioLayersP3", nBioLayersP3)
       call MPAS_pool_get_dimension(mesh, "TWO", TWO)
       call MPAS_pool_get_dimension(mesh, "nAlgae", nAlgae)
       call MPAS_pool_get_dimension(mesh, "nDOC", nDOC)
       call MPAS_pool_get_dimension(mesh, "nDIC", nDIC)
       call MPAS_pool_get_dimension(mesh, "nDON", nDON)
       call MPAS_pool_get_dimension(mesh, "nParticulateIron", nParticulateIron)
       call MPAS_pool_get_dimension(mesh, "nDissolvedIron", nDissolvedIron)

       call MPAS_pool_get_array(tracers, "verticalAerosolsConc",verticalAerosolsConc,1)
       call MPAS_pool_get_array(tracers, "verticalAerosolsSnow",verticalAerosolsSnow,1)
       call MPAS_pool_get_array(tracers, "verticalAerosolsIce",verticalAerosolsIce,1)
       call MPAS_pool_get_array(tracers, "verticalAlgaeConc",verticalAlgaeConc,1)
       call MPAS_pool_get_array(tracers, "verticalDOCConc",verticalDOCConc,1)
       call MPAS_pool_get_array(tracers, "verticalDICConc",verticalDICConc,1)
       call MPAS_pool_get_array(tracers, "verticalDONConc",verticalDONConc,1)
       call MPAS_pool_get_array(tracers, "verticalNitrateConc",verticalNitrateConc,1)
       call MPAS_pool_get_array(tracers, "verticalSilicateConc",verticalSilicateConc,1)
       call MPAS_pool_get_array(tracers, "verticalAmmoniumConc",verticalAmmoniumConc,1)
       call MPAS_pool_get_array(tracers, "verticalDMSConc",verticalDMSConc,1)
       call MPAS_pool_get_array(tracers, "verticalDMSPpConc",verticalDMSPpConc,1)
       call MPAS_pool_get_array(tracers, "verticalDMSPdConc",verticalDMSPdConc,1)
       call MPAS_pool_get_array(tracers, "verticalNonreactiveConc",verticalNonreactiveConc,1)
       call MPAS_pool_get_array(tracers, "verticalHumicsConc",verticalHumicsConc,1)
       call MPAS_pool_get_array(tracers, "verticalParticulateIronConc",verticalParticulateIronConc,1)
       call MPAS_pool_get_array(tracers, "verticalDissolvedIronConc",verticalDissolvedIronConc,1)
       call MPAS_pool_get_array(tracers, "verticalAlgaeSnow",verticalAlgaeSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDOCSnow",verticalDOCSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDICSnow",verticalDICSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDONSnow",verticalDONSnow,1)
       call MPAS_pool_get_array(tracers, "verticalNitrateSnow",verticalNitrateSnow,1)
       call MPAS_pool_get_array(tracers, "verticalSilicateSnow",verticalSilicateSnow,1)
       call MPAS_pool_get_array(tracers, "verticalAmmoniumSnow",verticalAmmoniumSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDMSSnow",verticalDMSSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDMSPpSnow",verticalDMSPpSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDMSPdSnow",verticalDMSPdSnow,1)
       call MPAS_pool_get_array(tracers, "verticalNonreactiveSnow",verticalNonreactiveSnow,1)
       call MPAS_pool_get_array(tracers, "verticalHumicsSnow",verticalHumicsSnow,1)
       call MPAS_pool_get_array(tracers, "verticalParticulateIronSnow",verticalParticulateIronSnow,1)
       call MPAS_pool_get_array(tracers, "verticalDissolvedIronSnow",verticalDissolvedIronSnow,1)
       call MPAS_pool_get_array(tracers, "verticalAlgaeIce",verticalAlgaeIce,1)
       call MPAS_pool_get_array(tracers, "verticalDOCIce",verticalDOCIce,1)
       call MPAS_pool_get_array(tracers, "verticalDICIce",verticalDICIce,1)
       call MPAS_pool_get_array(tracers, "verticalDONIce",verticalDONIce,1)
       call MPAS_pool_get_array(tracers, "verticalNitrateIce",verticalNitrateIce,1)
       call MPAS_pool_get_array(tracers, "verticalSilicateIce",verticalSilicateIce,1)
       call MPAS_pool_get_array(tracers, "verticalAmmoniumIce",verticalAmmoniumIce,1)
       call MPAS_pool_get_array(tracers, "verticalDMSIce",verticalDMSIce,1)
       call MPAS_pool_get_array(tracers, "verticalDMSPpIce",verticalDMSPpIce,1)
       call MPAS_pool_get_array(tracers, "verticalDMSPdIce",verticalDMSPdIce,1)
       call MPAS_pool_get_array(tracers, "verticalNonreactiveIce",verticalNonreactiveIce,1)
       call MPAS_pool_get_array(tracers, "verticalHumicsIce",verticalHumicsIce,1)
       call MPAS_pool_get_array(tracers, "verticalParticulateIronIce",verticalParticulateIronIce,1)
       call MPAS_pool_get_array(tracers, "verticalDissolvedIronIce",verticalDissolvedIronIce,1)

       do iCell = 1, nCellsSolve
          do iCategory = 1, nCategories

             ! aerosols
             if (config_use_zaerosols) then
                do iBioTracers = 1, nzAerosols
                   iSnowCount = (iBioTracers-1)*2
                   iIceCount = (iBioTracers-1)*nBioLayersP1
                   iBioData =  (iBioTracers-1)*nBioLayersP3

                   do iBioCount = 1,TWO
                      verticalAerosolsConc(iBioData+nBioLayersP1+iBioCount,iCategory,iCell) = &
                         verticalAerosolsSnow(iBioCount+iSnowCount,iCategory,iCell)
                   enddo
                   do iBioCount = 1, nBioLayersP1
                      verticalAerosolsConc(iBioData+iBioCount,iCategory,iCell) = &
                         verticalAerosolsIce(iBioCount+iIceCount,iCategory,iCell)
                   enddo
                enddo
             endif

             ! algal nitrogen
             if (config_use_vertical_biochemistry) then
                do iBioTracers = 1, nAlgae
                   iSnowCount = (iBioTracers-1)*2
                   iIceCount = (iBioTracers-1)*nBioLayersP1
                   iBioData =  (iBioTracers-1)*nBioLayersP3

                   do iBioCount = 1,TWO
                      verticalAlgaeConc(iBioData+nBioLayersP1+iBioCount,iCategory,iCell) = &
                         verticalAlgaeSnow(iBioCount+iSnowCount,iCategory,iCell)
                   enddo
                   do iBioCount = 1, nBioLayersP1
                      verticalAlgaeConc(iBioData+iBioCount,iCategory,iCell) = &
                         verticalAlgaeIce(iBioCount+iIceCount,iCategory,iCell)
                   enddo
                enddo
             endif

             ! dissolved organic and inorganic carbon
             if (config_use_carbon) then
                do iBioTracers = 1, nDOC
                   iSnowCount = (iBioTracers-1)*2
                   iIceCount = (iBioTracers-1)*nBioLayersP1
                   iBioData =  (iBioTracers-1)*nBioLayersP3

                   do iBioCount = 1,TWO
                      verticalDOCConc(iBioData+nBioLayersP1+iBioCount,iCategory,iCell) = &
                         verticalDOCSnow(iBioCount+iSnowCount,iCategory,iCell)
                   enddo
                   do iBioCount = 1, nBioLayersP1
                      verticalDOCConc(iBioData+iBioCount,iCategory,iCell) = &
                         verticalDOCIce(iBioCount+iIceCount,iCategory,iCell)
                   enddo
                enddo
                do iBioTracers = 1, nDIC
                   iSnowCount = (iBioTracers-1)*2
                   iIceCount = (iBioTracers-1)*nBioLayersP1
                   iBioData =  (iBioTracers-1)*nBioLayersP3

                   do iBioCount = 1,TWO
                      verticalDICConc(iBioData+nBioLayersP1+iBioCount,iCategory,iCell) = &
                         verticalDICSnow(iBioCount+iSnowCount,iCategory,iCell)
                   enddo
                   do iBioCount = 1, nBioLayersP1
                      verticalDICConc(iBioData+iBioCount,iCategory,iCell) = &
                         verticalDICIce(iBioCount+iIceCount,iCategory,iCell)
                   enddo
                enddo
             endif

             ! nitrate
             if (config_use_nitrate) then
                do iBioCount = 1,TWO
                   verticalNitrateConc(nBioLayersP1+iBioCount,iCategory,iCell) = &
                      verticalNitrateSnow(iBioCount,iCategory,iCell)
                enddo
                do iBioCount = 1, nBioLayersP1
                   verticalNitrateConc(iBioCount,iCategory,iCell) = &
                      verticalNitrateIce(iBioCount,iCategory,iCell)
                enddo
             endif

             ! ammonium
             if (config_use_ammonium) then
                do iBioCount = 1,TWO
                   verticalAmmoniumConc(nBioLayersP1+iBioCount,iCategory,iCell) = &
                      verticalAmmoniumSnow(iBioCount,iCategory,iCell)
                enddo
                do iBioCount = 1, nBioLayersP1
                   verticalAmmoniumConc(iBioCount,iCategory,iCell) = &
                      verticalAmmoniumIce(iBioCount,iCategory,iCell)
                enddo
             endif

             ! silicate
             if (config_use_silicate) then
                do iBioCount = 1,TWO
                   verticalSilicateConc(nBioLayersP1+iBioCount,iCategory,iCell) = &
                      verticalSilicateSnow(iBioCount,iCategory,iCell)
                enddo
                do iBioCount = 1, nBioLayersP1
                   verticalSilicateConc(iBioCount,iCategory,iCell) = &
                      verticalSilicateIce(iBioCount,iCategory,iCell)
                enddo
             endif

             ! DMS, DMSPp, DMPSd
             if (config_use_DMS) then
                do iBioCount = 1,TWO
                   verticalDMSConc(nBioLayersP1+iBioCount,iCategory,iCell) = &
                      verticalDMSSnow(iBioCount,iCategory,iCell)
                   verticalDMSPpConc(nBioLayersP1+iBioCount,iCategory,iCell) = &
                      verticalDMSPpSnow(iBioCount,iCategory,iCell)
                   verticalDMSPdConc(nBioLayersP1+iBioCount,iCategory,iCell) = &
                      verticalDMSPdSnow(iBioCount,iCategory,iCell)
                enddo
                do iBioCount = 1, nBioLayersP1
                   verticalDMSConc(iBioCount,iCategory,iCell) = &
                      verticalDMSIce(iBioCount,iCategory,iCell)
                   verticalDMSPpConc(iBioCount,iCategory,iCell) = &
                      verticalDMSPpIce(iBioCount,iCategory,iCell)
                   verticalDMSPdConc(iBioCount,iCategory,iCell) = &
                      verticalDMSPdIce(iBioCount,iCategory,iCell)
                enddo
             endif

             ! nonreactive tracer
             if (config_use_nonreactive) then
                do iBioCount = 1,TWO
                   verticalNonreactiveConc(nBioLayersP1+iBioCount,iCategory,iCell) = &
                      verticalNonreactiveSnow(iBioCount,iCategory,iCell)
                enddo
                do iBioCount = 1, nBioLayersP1
                   verticalNonreactiveConc(iBioCount,iCategory,iCell) = &
                      verticalNonreactiveIce(iBioCount,iCategory,iCell)
                enddo
             endif

             ! humics
             if (config_use_humics) then
                do iBioCount = 1,TWO
                   verticalHumicsConc(nBioLayersP1+iBioCount,iCategory,iCell) = &
                      verticalHumicsSnow(iBioCount,iCategory,iCell)
                enddo
                do iBioCount = 1, nBioLayersP1
                   verticalHumicsConc(iBioCount,iCategory,iCell) = &
                      verticalHumicsIce(iBioCount,iCategory,iCell)
                enddo
             endif

             ! proteins and amino acids
             if (config_use_DON) then
                do iBioTracers = 1, nDON
                   iSnowCount = (iBioTracers-1)*2
                   iIceCount = (iBioTracers-1)*nBioLayersP1
                   iBioData =  (iBioTracers-1)*nBioLayersP3

                   do iBioCount = 1,TWO
                      verticalDONConc(iBioData+nBioLayersP1+iBioCount,iCategory,iCell) = &
                         verticalDONSnow(iBioCount+iSnowCount,iCategory,iCell)
                   enddo
                   do iBioCount = 1, nBioLayersP1
                      verticalDONConc(iBioData+iBioCount,iCategory,iCell) = &
                         verticalDONIce(iBioCount+iIceCount,iCategory,iCell)
                   enddo
                enddo
             endif

             ! particulate and dissolved iron
             if (config_use_iron) then
                do iBioTracers = 1, nParticulateIron
                   iSnowCount = (iBioTracers-1)*2
                   iIceCount = (iBioTracers-1)*nBioLayersP1
                   iBioData =  (iBioTracers-1)*nBioLayersP3

                   do iBioCount = 1,TWO
                      verticalDissolvedIronConc(iBioData+nBioLayersP1+iBioCount,iCategory,iCell) = &
                         verticalDissolvedIronSnow(iBioCount+iSnowCount,iCategory,iCell)
                   enddo
                   do iBioCount = 1, nBioLayersP1
                      verticalDissolvedIronConc(iBioData+iBioCount,iCategory,iCell) = &
                         verticalDissolvedIronIce(iBioCount+iIceCount,iCategory,iCell)
                   enddo
                enddo
             endif

          enddo
       enddo

       block => block % next
    end do

  end subroutine column_combine_snow_ice_tracers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_ocean_carbon_flux
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 26 May 2020
!> \details Calculate the ocean carbon flux
!> by summing the appropriate biogeochemical tracer fluxes in units of mmol C/m2/s
!>
!>      ocean carbon flux = algal nitrogen group fluxes * (C to N ratios)
!>                   + dissolved carbon group fluxes
!>                   + dissolved organic nitrogen * (C to N ratio)
!>                   + dissolved inorganic carbon fluxes + humic fluxes
!
!-----------------------------------------------------------------------

  subroutine seaice_ocean_carbon_flux(block,oceanCarbonFlux,oceanBioFluxes,iCell)

    real(kind=RKIND), dimension(:), intent(out) :: &
         oceanCarbonFlux

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         oceanBioFluxes

    integer, intent(in) :: &
         iCell

    type(block_type), intent(in) :: &
         block

    logical, pointer :: &
         config_use_column_biogeochemistry, &
         config_use_vertical_biochemistry, &
         config_use_carbon, &
         config_use_DON, &
         config_use_humics

    integer, pointer :: &
         nAlgae, &
         nDOC, &
         nDON, &
         nDIC

    type(MPAS_pool_type), pointer :: &
         mesh, &
         biogeochemistry

    real(kind=RKIND), pointer :: &
         config_ratio_C_to_N_diatoms, &
         config_ratio_C_to_N_small_plankton, &
         config_ratio_C_to_N_phaeocystis, &
         config_ratio_C_to_N_proteins

    integer, pointer :: &
         nCategories, &
         nZBGCTracers, &
         maxAlgaeType, &
         maxDOCType, &
         maxDICType, &
         maxDONType, &
         maxIronType, &
         maxBCType, &
         maxDustType, &
         maxAerosolType

    real(kind=RKIND), dimension(:), allocatable :: &
         ratio_C_to_N

    real(kind=RKIND), dimension(:,:), allocatable :: &
         oceanBioFluxesAll

    integer :: &
         iBioTracers, &
         iBioData, &
         iCategory

    call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(block % configs, "config_use_DON", config_use_DON)
    call MPAS_pool_get_config(block % configs, "config_use_humics",config_use_humics)
    call MPAS_pool_get_config(block % configs, "config_ratio_C_to_N_diatoms", config_ratio_C_to_N_diatoms)
    call MPAS_pool_get_config(block % configs, "config_ratio_C_to_N_small_plankton", config_ratio_C_to_N_small_plankton)
    call MPAS_pool_get_config(block % configs, "config_ratio_C_to_N_phaeocystis", config_ratio_C_to_N_phaeocystis)
    call MPAS_pool_get_config(block % configs, "config_ratio_C_to_N_proteins", config_ratio_C_to_N_proteins)

    call MPAS_pool_get_dimension(block % dimensions, "nAlgae", nAlgae)
    call MPAS_pool_get_dimension(block % dimensions, "nDOC", nDOC)
    call MPAS_pool_get_dimension(block % dimensions, "nDIC", nDIC)
    call MPAS_pool_get_dimension(block % dimensions, "nDON", nDON)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)

    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
    call MPAS_pool_get_dimension(mesh, "nZBGCTracers", nZBGCTracers)
    call MPAS_pool_get_dimension(mesh, "maxAlgaeType", maxAlgaeType)
    call MPAS_pool_get_dimension(mesh, "maxDOCType", maxDOCType)
    call MPAS_pool_get_dimension(mesh, "maxDICType", maxDICType)
    call MPAS_pool_get_dimension(mesh, "maxDONType", maxDONType)
    call MPAS_pool_get_dimension(mesh, "maxAerosolType", maxAerosolType)
    call MPAS_pool_get_dimension(mesh, "maxIronType", maxIronType)
    call MPAS_pool_get_dimension(mesh, "maxBCType", maxBCType)
    call MPAS_pool_get_dimension(mesh, "maxDustType", maxDustType)


    allocate(oceanBioFluxesAll(nZBGCTracers,nCategories))
    allocate(ratio_C_to_N(3))

    ratio_C_to_N(1) = config_ratio_C_to_N_diatoms
    ratio_C_to_N(2) = config_ratio_C_to_N_small_plankton
    ratio_C_to_N(3) = config_ratio_C_to_N_phaeocystis

    if (config_use_column_biogeochemistry) then

       do iCategory = 1, nCategories

       oceanCarbonFlux(iCategory)   = 0.0_RKIND
       oceanBioFluxesAll(:,iCategory) = 0.0_RKIND

       do iBioTracers = 1, ciceTracerObject % nBioTracers
          iBioData = ciceTracerObject % index_LayerIndexToDataArray(iBioTracers)
          oceanBioFluxesAll(iBioData,iCategory) = oceanBioFluxes(iBioTracers,iCategory,iCell)
       enddo
       iBioData = 0

       ! Algae
       do iBioTracers = 1, maxAlgaeType
          iBioData = iBioData+1
          oceanCarbonFlux(iCategory) = oceanCarbonFlux(iCategory) + &
             oceanBioFluxesAll(iBioData,iCategory) * ratio_C_to_N(iBioTracers)
       enddo

       ! Nitrate
       iBioData = iBioData+1

       ! Polysaccharids and Lipids
       do iBioTracers = 1, maxDOCType
          iBioData = iBioData+1
          oceanCarbonFlux(iCategory) = oceanCarbonFlux(iCategory) + &
             oceanBioFluxesAll(iBioData,iCategory)
       enddo

       ! DIC
       do iBioTracers = 1, maxDICType
          iBioData = iBioData+1
          oceanCarbonFlux(iCategory) = oceanCarbonFlux(iCategory) + &
             oceanBioFluxesAll(iBioData,iCategory)
       enddo

       ! + Chlorophyll (maxAlgaeType) + Ammonium (1) + Silicate (1) + DMSPp (1) + DMSPd (1)
       ! + DMS (1) + PON (1)

       iBioData = iBioData+maxAlgaeType + 6

       ! DON
       do iBioTracers = 1, maxDONType
          iBioData = iBioData+1
          oceanCarbonFlux(iCategory) = oceanCarbonFlux(iCategory) + &
             oceanBioFluxesAll(iBioData,iCategory) * config_ratio_C_to_N_proteins
       enddo

       ! + dFe (maxIronType) + pFe (maxIronType)
       ! + Black Carbon (maxBCType) + Dust (maxDustType)

       iBioData = iBioData + 2*maxIronType + maxBCType + maxDustType

       ! Humics
       iBioData = iBioData+1
       oceanCarbonFlux(iCategory) = oceanCarbonFlux(iCategory) + &
          oceanBioFluxesAll(iBioData,iCategory)

       enddo ! nCategories
    endif

    deallocate(oceanBioFluxesAll)
    deallocate(ratio_C_to_N)

  end subroutine seaice_ocean_carbon_flux


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_ocean_carbon_flux_cell
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 26 May 2020
!> \details Calculate the ocean carbon flux
!> by summing the appropriate biogeochemical tracer fluxes in units of mmol C/m2/s
!>
!>      ocean carbon flux = algal nitrogen group fluxes * (C to N ratios)
!>                   + dissolved carbon group fluxes
!>                   + dissolved organic nitrogen * (C to N ratio)
!>                   + dissolved inorganic carbon fluxes + humic fluxes
!
!-----------------------------------------------------------------------

  subroutine seaice_ocean_carbon_flux_cell(block,oceanCarbonFlux,oceanBioFluxes)

    real(kind=RKIND), intent(out) :: &
         oceanCarbonFlux

    real(kind=RKIND), dimension(:), intent(in) :: &
         oceanBioFluxes

    type(block_type), intent(in) :: &
         block

    logical, pointer :: &
         config_use_column_biogeochemistry, &
         config_use_vertical_biochemistry, &
         config_use_carbon, &
         config_use_DON, &
         config_use_humics

    integer, pointer :: &
         nAlgae, &
         nDOC, &
         nDON, &
         nDIC

    type(MPAS_pool_type), pointer :: &
         mesh, &
         biogeochemistry

    real(kind=RKIND), pointer :: &
         config_ratio_C_to_N_diatoms, &
         config_ratio_C_to_N_small_plankton, &
         config_ratio_C_to_N_phaeocystis, &
         config_ratio_C_to_N_proteins

    integer, pointer :: &
         nZBGCTracers, &
         maxAlgaeType, &
         maxDOCType, &
         maxDICType, &
         maxDONType, &
         maxIronType, &
         maxBCType, &
         maxDustType, &
         maxAerosolType

    real(kind=RKIND), dimension(:), allocatable :: &
         ratio_C_to_N

    real(kind=RKIND), dimension(:), allocatable :: &
         oceanBioFluxesAll

    integer :: &
         iBioTracers, &
         iBioData

    call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(block % configs, "config_use_DON", config_use_DON)
    call MPAS_pool_get_config(block % configs, "config_use_humics",config_use_humics)
    call MPAS_pool_get_config(block % configs, "config_ratio_C_to_N_diatoms", config_ratio_C_to_N_diatoms)
    call MPAS_pool_get_config(block % configs, "config_ratio_C_to_N_small_plankton", config_ratio_C_to_N_small_plankton)
    call MPAS_pool_get_config(block % configs, "config_ratio_C_to_N_phaeocystis", config_ratio_C_to_N_phaeocystis)
    call MPAS_pool_get_config(block % configs, "config_ratio_C_to_N_proteins", config_ratio_C_to_N_proteins)

    call MPAS_pool_get_dimension(block % dimensions, "nAlgae", nAlgae)
    call MPAS_pool_get_dimension(block % dimensions, "nDOC", nDOC)
    call MPAS_pool_get_dimension(block % dimensions, "nDIC", nDIC)
    call MPAS_pool_get_dimension(block % dimensions, "nDON", nDON)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)

    call MPAS_pool_get_dimension(mesh, "nZBGCTracers", nZBGCTracers)
    call MPAS_pool_get_dimension(mesh, "maxAlgaeType", maxAlgaeType)
    call MPAS_pool_get_dimension(mesh, "maxDOCType", maxDOCType)
    call MPAS_pool_get_dimension(mesh, "maxDICType", maxDICType)
    call MPAS_pool_get_dimension(mesh, "maxDONType", maxDONType)
    call MPAS_pool_get_dimension(mesh, "maxAerosolType", maxAerosolType)
    call MPAS_pool_get_dimension(mesh, "maxIronType", maxIronType)
    call MPAS_pool_get_dimension(mesh, "maxBCType", maxBCType)
    call MPAS_pool_get_dimension(mesh, "maxDustType", maxDustType)

    allocate(oceanBioFluxesAll(nZBGCTracers))
    allocate(ratio_C_to_N(3))

    ratio_C_to_N(1) = config_ratio_C_to_N_diatoms
    ratio_C_to_N(2) = config_ratio_C_to_N_small_plankton
    ratio_C_to_N(3) = config_ratio_C_to_N_phaeocystis

    if (config_use_column_biogeochemistry) then

       oceanCarbonFlux   = 0.0_RKIND
       oceanBioFluxesAll(:) = 0.0_RKIND

       do iBioTracers = 1, ciceTracerObject % nBioTracers
          iBioData = ciceTracerObject % index_LayerIndexToDataArray(iBioTracers)
          oceanBioFluxesAll(iBioData) = oceanBioFluxes(iBioTracers)
       enddo
       iBioData = 0

       ! Algae
       do iBioTracers = 1, maxAlgaeType
          iBioData = iBioData+1
          oceanCarbonFlux = oceanCarbonFlux + &
             oceanBioFluxesAll(iBioData) * ratio_C_to_N(iBioTracers)
       enddo

       ! Nitrate
       iBioData = iBioData+1

       ! Polysaccharids and Lipids
       do iBioTracers = 1, maxDOCType
          iBioData = iBioData+1
          oceanCarbonFlux = oceanCarbonFlux + &
             oceanBioFluxesAll(iBioData)
       enddo

       ! DIC
       do iBioTracers = 1, maxDICType
          iBioData = iBioData+1
          oceanCarbonFlux = oceanCarbonFlux + &
             oceanBioFluxesAll(iBioData)
       enddo

       ! + Chlorophyll (maxAlgaeType) + Ammonium (1) + Silicate (1) + DMSPp (1) + DMSPd (1)
       ! + DMS (1) + PON (1)

       iBioData = iBioData+maxAlgaeType + 6

       ! DON
       do iBioTracers = 1, maxDONType
          iBioData = iBioData+1
          oceanCarbonFlux = oceanCarbonFlux + &
             oceanBioFluxesAll(iBioData) * config_ratio_C_to_N_proteins
       enddo

       ! + dFe (maxIronType) + pFe (maxIronType)
       ! + Black Carbon (maxBCType) + Dust (maxDustType)

       iBioData = iBioData + 2*maxIronType + maxBCType + maxDustType

       ! Humics
       iBioData = iBioData+1
       oceanCarbonFlux = oceanCarbonFlux + &
          oceanBioFluxesAll(iBioData)

    endif

    deallocate(oceanBioFluxesAll)
    deallocate(ratio_C_to_N)

  end subroutine seaice_ocean_carbon_flux_cell


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_total_carbon_content_category
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 26 May 2020
!> \details Calculate the total carbon concentration in the sea ice category
!> by summing the appropriate biogeochemical tracers in units of mmol C
!>
!>      Total carbon = algal nitrogen groups * (C to N ratios) + dissolved carbon groups
!>                   + dissolved inorganic carbon + humic material
!>                   + dissolved organic nitrogen * (C to N ratio)
!
!-----------------------------------------------------------------------

  subroutine seaice_total_carbon_content_category(block,totalCarbonContentCategory,iceAreaCategory,iceVolumeCategory,iCell)

    use seaice_constants, only: &
         skeletalLayerThickness, &
         seaicePuny

    real(kind=RKIND), dimension(:), intent(out) :: &
         totalCarbonContentCategory

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         iceAreaCategory, &
         iceVolumeCategory

    integer, intent(in) :: &
         iCell

    type(block_type), intent(in) :: &
         block

    logical, pointer :: &
         config_use_skeletal_biochemistry, &
         config_use_vertical_biochemistry, &
         config_use_vertical_tracers, &
         config_use_carbon, &
         config_use_DON, &
         config_use_humics

    integer, pointer :: &
         nCategories, &
         nBioLayersP1, &
         nBioLayers, &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON

    type(MPAS_pool_type), pointer :: &
         mesh, &
         biogeochemistry, &
         tracers

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         skeletalAlgaeConc, &
         skeletalDOCConc, &
         skeletalDICConc, &
         skeletalDONConc, &
         skeletalHumicsConc, &
         verticalAlgaeConc, &
         verticalDOCConc, &
         verticalDICConc, &
         verticalDONConc, &
         verticalHumicsConc, &
         brineFraction

    real(kind=RKIND), pointer :: &
         config_ratio_C_to_N_diatoms, &
         config_ratio_C_to_N_small_plankton, &
         config_ratio_C_to_N_phaeocystis, &
         config_ratio_C_to_N_proteins

    real(kind=RKIND), dimension(:), allocatable :: &
         ratio_C_to_N, &
         verticalGridSpace

    real(kind=RKIND) :: &
         brineHeight

    integer :: &
         iBioTracers, &
         iBioCount, &
         iLayers, &
         iCategory

    call MPAS_pool_get_config(block % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_vertical_tracers", config_use_vertical_tracers)
    call MPAS_pool_get_config(block % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(block % configs, "config_use_DON", config_use_DON)
    call MPAS_pool_get_config(block % configs, "config_use_humics",config_use_humics)
    call MPAS_pool_get_config(block % configs, "config_ratio_C_to_N_diatoms", config_ratio_C_to_N_diatoms)
    call MPAS_pool_get_config(block % configs, "config_ratio_C_to_N_small_plankton", config_ratio_C_to_N_small_plankton)
    call MPAS_pool_get_config(block % configs, "config_ratio_C_to_N_phaeocystis", config_ratio_C_to_N_phaeocystis)
    call MPAS_pool_get_config(block % configs, "config_ratio_C_to_N_proteins", config_ratio_C_to_N_proteins)

    call MPAS_pool_get_dimension(block % dimensions, "nBioLayers", nBioLayers)
    call MPAS_pool_get_dimension(block % dimensions, "nBioLayersP1", nBioLayersP1)
    call MPAS_pool_get_dimension(block % dimensions, "nAlgae", nAlgae)
    call MPAS_pool_get_dimension(block % dimensions, "nDOC", nDOC)
    call MPAS_pool_get_dimension(block % dimensions, "nDIC", nDIC)
    call MPAS_pool_get_dimension(block % dimensions, "nDON", nDON)

    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)
    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(tracers, "skeletalAlgaeConc", skeletalAlgaeConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDOCConc", skeletalDOCConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDICConc", skeletalDICConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalDONConc", skeletalDONConc, 1)
    call MPAS_pool_get_array(tracers, "skeletalHumicsConc", skeletalHumicsConc, 1)
    call MPAS_pool_get_array(tracers, "verticalAlgaeConc", verticalAlgaeConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDOCConc", verticalDOCConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDICConc", verticalDICConc, 1)
    call MPAS_pool_get_array(tracers, "verticalDONConc", verticalDONConc, 1)
    call MPAS_pool_get_array(tracers, "verticalHumicsConc", verticalHumicsConc, 1)
    call MPAS_pool_get_array(tracers, "brineFraction", brineFraction, 1)

    allocate(ratio_C_to_N(3))
    allocate(verticalGridSpace(nBioLayersP1))

    ratio_C_to_N(1) = config_ratio_C_to_N_diatoms
    ratio_C_to_N(2) = config_ratio_C_to_N_small_plankton
    ratio_C_to_N(3) = config_ratio_C_to_N_phaeocystis


    verticalGridSpace(:) = 1.0_RKIND/real(nBioLayers,kind=RKIND)
    verticalGridSpace(1) = verticalGridSpace(1)/2.0_RKIND
    verticalGridSpace(nBioLayersP1) = verticalGridSpace(1)
    totalCarbonContentCategory(:) = 0.0_RKIND


    if (config_use_skeletal_biochemistry) then

       do iCategory = 1, nCategories
       ! algal nitrogen
       do iBioTracers = 1, nAlgae
          totalCarbonContentCategory(iCategory) = totalCarbonContentCategory(iCategory) +  skeletalAlgaeConc(iBioTracers,iCategory,iCell)* &
             skeletalLayerThickness * ratio_C_to_N(iBioTracers)
       enddo

       if (config_use_carbon) then
          ! DOC
          do iBioTracers = 1, nDOC
             totalCarbonContentCategory(iCategory) = totalCarbonContentCategory(iCategory) +  skeletalDOCConc(iBioTracers,iCategory,iCell)* &
                skeletalLayerThickness
          enddo

          ! DIC
          do iBioTracers = 1, nDIC
             totalCarbonContentCategory(iCategory) = totalCarbonContentCategory(iCategory) +  skeletalDICConc(iBioTracers,iCategory,iCell)* &
                skeletalLayerThickness
          enddo
       endif

       if (config_use_DON) then
          ! DON
          do iBioTracers = 1, nDON
             totalCarbonContentCategory(iCategory) = totalCarbonContentCategory(iCategory) +  skeletalDONConc(iBioTracers,iCategory,iCell)* &
                config_ratio_C_to_N_proteins * skeletalLayerThickness
          enddo
       endif

       ! humic material
       if (config_use_humics) &
          totalCarbonContentCategory(iCategory) = totalCarbonContentCategory(iCategory) +  skeletalHumicsConc(1,iCategory,iCell)* &
             skeletalLayerThickness
        enddo
    elseif (config_use_vertical_tracers) then

       do iCategory = 1, nCategories
          brineHeight = 0.0_RKIND
          if (iceAreaCategory(iCategory,iCell) > seaicePuny) then
             brineHeight = iceVolumeCategory(iCategory,iCell)/iceAreaCategory(iCategory,iCell) * brineFraction(1,iCategory,iCell)
          endif

          if (config_use_vertical_biochemistry) then
             iBioCount = 0
             ! algal nitrogen
             do iBioTracers = 1, nAlgae
                do iLayers = 1,nBioLayersP1
                   iBiocount = iBiocount + 1
                   totalCarbonContentCategory(iCategory) = totalCarbonContentCategory(iCategory) + &
                      verticalAlgaeConc(iBioCount,iCategory,iCell) * ratio_C_to_N(iBioTracers) * &
                      verticalGridSpace(iLayers) * brineHeight
                 enddo
                 iBioCount = iBioCount+2   ! snow layers
              enddo
           endif

           if (config_use_carbon) then
              iBioCount = 0
              ! DOC
              do iBioTracers = 1, nDOC
                 do iLayers = 1,nBioLayersP1
                    iBioCount = iBioCount + 1
                    totalCarbonContentCategory(iCategory) = totalCarbonContentCategory(iCategory) + &
                       verticalDOCConc(iBioCount,iCategory,iCell) * verticalGridSpace(iLayers) * brineHeight
                 enddo
                 iBioCount = iBioCount+2  ! snow layers
              enddo
              iBioCount = 0
              ! DIC
              do iBioTracers = 1, nDIC

                 do iLayers = 1,nBioLayersP1
                    iBioCount = iBioCount + 1
                    totalCarbonContentCategory(iCategory) = totalCarbonContentCategory(iCategory) + &
                       verticalDICConc(iBioCount,iCategory,iCell)  * verticalGridSpace(iLayers) * brineHeight
                 enddo
                 iBioCount = iBioCount + 2 ! snow layers
              enddo
           endif

          if (config_use_DON) then
             iBioCount = 0
             ! dissolved organic nitrogen
             do iBioTracers = 1, nDON
                do iLayers = 1,nBioLayersP1
                   iBiocount = iBiocount + 1
                   totalCarbonContentCategory(iCategory) = totalCarbonContentCategory(iCategory) + &
                      verticalDONConc(iBioCount,iCategory,iCell) * config_ratio_C_to_N_proteins * &
                      verticalGridSpace(iLayers) * brineHeight
                 enddo
                 iBioCount = iBioCount+2   ! snow layers
              enddo
           endif

           ! humic material
           if (config_use_humics) then
              do iLayers = 1, nBioLayersP1
                 totalCarbonContentCategory(iCategory) = totalCarbonContentCategory(iCategory) + &
                    verticalHumicsConc(iLayers,iCategory,iCell)  * verticalGridSpace(iLayers) * brineHeight
              enddo
           endif
       enddo
    endif

    deallocate(ratio_C_to_N)
    deallocate(verticalGridSpace)

  end subroutine seaice_total_carbon_content_category

!-----------------------------------------------------------------------

  subroutine seaice_icepack_write_warnings(logAsErrors)

    use icepack_intfc, only: &
         icepack_warnings_clear, &
         icepack_warnings_getall

    character(len=256), dimension(:), allocatable :: & ! icepack char_len_long=256
         warnings

    logical, intent(in) :: &
         logAsErrors

    integer :: &
         iWarning

    call icepack_warnings_getall(warnings)

    if (logAsErrors) then
       do iWarning = 1, size(warnings)
          call mpas_log_write(trim(warnings(iWarning)), messageType=MPAS_LOG_ERR)
       enddo                    ! iWarning
       call mpas_log_write("icepack aborted", MPAS_LOG_CRIT)
    else
       do iWarning = 1, size(warnings)
          call mpas_log_write(trim(warnings(iWarning)), messageType=MPAS_LOG_WARN)
       enddo                    ! iWarning
    endif

    call icepack_warnings_clear()

  end subroutine seaice_icepack_write_warnings

  !-----------------------------------------------------------------------

end module seaice_icepack
