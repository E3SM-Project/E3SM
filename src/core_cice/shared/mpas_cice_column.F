#define COMMA ,
#define COLUMN_DEBUG_WRITE(M) ! write(stderrUnit,*) M
#define COLUMN_WARNING_WRITE(M) write(stderrUnit,*) 'WARNING: '//M
#define COLUMN_ERROR_WRITE(M) write(stderrUnit,*) 'ERROR: '//M

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

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_timer

  use cice_debug

  use ice_kinds_mod, only: &
       char_len_long

  implicit none

  private
  save

  public :: &
       cice_init_column_physics_package_parameters, &
       cice_init_column_physics_package_variables, &
       cice_column_predynamics_time_integration, &
       cice_column_dynamics_time_integration, &
       cice_column_postdynamics_time_integration, &
       cice_init_column_shortwave, &
       cice_column_aggregate, &
       cice_column_initial_air_drag_coefficient, &
       cice_column_reinitialize_fluxes, &
       cice_column_reinitialize_diagnostics_thermodynamics, &
       cice_column_reinitialize_diagnostics_dynamics, &
       cice_column_reinitialize_diagnostics_bgc, &
       cice_column_coupling_prep, &
       cice_column_finalize

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

     ! category tracer array
     real(kind=RKIND), dimension(:,:), allocatable :: tracerArrayCategory ! trcrn

     ! cell tracer array
     real(kind=RKIND), dimension(:), allocatable :: tracerArrayCell ! trcr

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
          index_aerosols                 ! nt_aero

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

     ! indexes of bio tracers with types in tracer array
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

  ! warnings string kind
  integer, parameter :: strKINDWarnings = char_len_long

contains

!-----------------------------------------------------------------------
! Initialize Column Physics Package
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_column_physics_package_parameters
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 18th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_init_column_physics_package_parameters(domain)

    type(domain_type), intent(inout) :: domain

    logical, pointer :: &
         config_use_column_package

    call MPAS_pool_get_config(domain % configs, "config_use_column_package", config_use_column_package)
    if (config_use_column_package) then

       ! set non activated variable pointers to other memory
       call init_column_non_activated_pointers(domain)

       ! initialize the column package tracer object
       call init_column_tracer_object(domain, ciceTracerObject)

       ! initialize the column package parameters
       call init_column_package_parameters(domain, ciceTracerObject)

    endif

  end subroutine cice_init_column_physics_package_parameters

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_column_physics_package_variables
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 18th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_init_column_physics_package_variables(domain, clock)

    type(domain_type), intent(inout) :: domain

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    logical, pointer :: &
         config_use_column_package, &
         config_do_restart, &
         config_use_column_biogeochemistry, &
         config_use_column_shortwave

    call MPAS_pool_get_config(domain % configs, "config_use_column_package", config_use_column_package)
    call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

    if (config_use_column_package) then

       ! initialize level ice tracers
       call init_column_level_ice_tracers(domain)

       ! initialize the itd thickness classes
       call init_column_itd(domain)

       ! initialize thermodynamic tracer profiles
       call init_column_thermodynamic_profiles(domain)

       ! initialize biogoechemistry profiles
       if (config_use_column_biogeochemistry) &
            call init_column_biogeochemistry_profiles(domain, ciceTracerObject)

       ! history variables
       call init_column_history_variables(domain)

       ! shortwave
       call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)
       call MPAS_pool_get_config(domain % configs, "config_use_column_shortwave", config_use_column_shortwave)
       if (config_do_restart .and. config_use_column_shortwave) &
            call cice_init_column_shortwave(domain, clock)

    endif

  end subroutine cice_init_column_physics_package_variables

!-----------------------------------------------------------------------
! column package initialization routine wrappers
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_itd
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 5th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_itd(domain)

    use ice_colpkg, only: colpkg_init_itd

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         initial

    real(kind=RKIND), dimension(:), pointer :: &
         categoryThicknessLimits

    integer, pointer :: &
         nCategories

    logical :: &
         abortFlag

    character(len=strKIND) :: &
         abortMessage

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nCategories", nCategories)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "initial", initial)

       call MPAS_pool_get_array(initial, "categoryThicknessLimits", categoryThicknessLimits)

       call colpkg_init_itd(&
            nCategories, &
            categoryThicknessLimits, &
            abortFlag, &
            abortMessage)

       ! code abort
       if (abortFlag) then
          COLUMN_ERROR_WRITE("init_column_itd")
          COLUMN_ERROR_WRITE(trim(abortMessage))
          call MPAS_dmpar_global_abort("MPAS-seaice: init_column_itd: "//trim(abortMessage))
       endif

       block => block % next
    end do

  end subroutine init_column_itd

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_thermo
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 5th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_thermodynamic_profiles(domain)

    use ice_colpkg, only: &
         colpkg_init_thermo, &
         colpkg_liquidus_temperature

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

       call colpkg_init_thermo(&
            nIceLayers, &
            initialSalinityProfileVertical)

       call MPAS_pool_get_subpool(block % structs, "initial", initial)

       call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)
       call MPAS_pool_get_array(initial, "initialMeltingTemperatureProfile", initialMeltingTemperatureProfile)

       do iCell = 1, nCellsSolve
          do iIceLayer = 1, nIceLayers + 1

             ! these profiles are not used by mushy
             initialSalinityProfile(iIceLayer,iCell)           = initialSalinityProfileVertical(iIceLayer)
             initialMeltingTemperatureProfile(iIceLayer,iCell) = &
                  colpkg_liquidus_temperature(initialSalinityProfileVertical(iIceLayer))

          enddo ! iIceLayer
       enddo ! iCell

       deallocate(initialSalinityProfileVertical)

       block => block % next
    end do

  end subroutine init_column_thermodynamic_profiles

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_shortwave
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 5th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_init_column_shortwave(domain, clock)

    use ice_colpkg, only: &
         colpkg_init_orbit, &
         colpkg_get_warnings, &
         colpkg_clear_warnings

    use cice_constants, only: &
         cicePuny

    type(domain_type), intent(inout) :: domain

    type(MPAS_clock_type), intent(in) :: clock

    logical :: &
         abortFlag

    character(len=strKIND) :: &
         abortMessage

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

    character(len=strKINDWarnings), dimension(:), allocatable :: &
         warnings

    call MPAS_pool_get_config(domain % configs, "config_shortwave_type", config_shortwave_type)
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    if (trim(config_shortwave_type) == "dEdd") then

       abortFlag = .false.
       abortMessage = ""

       call colpkg_clear_warnings()
       call colpkg_init_orbit(&
            abortFlag, &
            abortMessage)
       call colpkg_get_warnings(warnings)
       call column_write_warnings(warnings)

       if (abortFlag) then
          COLUMN_ERROR_WRITE("colpkg_init_orbit")
          COLUMN_ERROR_WRITE(trim(abortMessage))
          call MPAS_dmpar_global_abort("MPAS-seaice: "//trim(abortMessage))
      endif

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
             if (iceAreaCategory(1,iCategory,iCell) > cicePuny) then

                albedoVisibleDirectCell(iCell)  = albedoVisibleDirectCell(iCell)  + &
                     albedoVisibleDirectCategory(iCategory,iCell)  * iceAreaCategory(1,iCategory,iCell)
                albedoVisibleDiffuseCell(iCell) = albedoVisibleDiffuseCell(iCell) + &
                     albedoVisibleDiffuseCategory(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)
                albedoIRDirectCell(iCell)       = albedoIRDirectCell(iCell)       + &
                     albedoIRDirectCategory(iCategory,iCell)       * iceAreaCategory(1,iCategory,iCell)
                albedoIRDiffuseCell(iCell)      = albedoIRDiffuseCell(iCell)      + &
                     albedoIRDiffuseCategory(iCategory,iCell)      * iceAreaCategory(1,iCategory,iCell)

                if (solarZenithAngleCosine(iCell) > cicePuny) then ! sun above horizon

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

  end subroutine cice_init_column_shortwave

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_thermodynamic_tracers
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 5th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_thermodynamic_tracers(domain)

    use ice_colpkg, only: colpkg_init_trcr

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
         nIceLayers, &
         nSnowLayers, &
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
       call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)

       call MPAS_pool_get_array(atmos_coupling, "airTemperature", airTemperature)

       call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)

       call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)
       call MPAS_pool_get_array(initial, "initialMeltingTemperatureProfile", initialMeltingTemperatureProfile)

       call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
       call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy, 1)
       call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy, 1)

       do iCell = 1, nCellsSolve
          do iCategory = 1, nCategories

             call colpkg_init_trcr(&
                  airTemperature(iCell), &
                  seaFreezingTemperature(iCell), &
                  initialSalinityProfile(:,iCell), &
                  initialMeltingTemperatureProfile(:,iCell), &
                  surfaceTemperature(1,iCategory,iCell), &
                  nIceLayers, &
                  nSnowLayers, &
                  iceEnthalpy(:,iCategory,iCell), &
                  snowEnthalpy(:,iCategory,iCell))

          enddo ! iCategory
       enddo ! iCell

       block => block % next
    end do

  end subroutine init_column_thermodynamic_tracers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_level_ice_tracers
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 6th Feburary 2015
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
!  cice_column_finalize
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 29th October 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_column_finalize(domain)

    type(domain_type), intent(inout) :: domain

    call finalize_column_non_activated_pointers(domain)

  end subroutine cice_column_finalize

!-----------------------------------------------------------------------
! runtime
!-----------------------------------------------------------------------
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_predynamics_time_integration
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 6th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_column_predynamics_time_integration(domain, clock)

    type(domain_type), intent(inout) :: domain

    type(MPAS_clock_type), intent(in) :: clock

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_shortwave, &
         config_use_column_vertical_thermodynamics, &
         config_use_column_biogeochemistry, &
         config_use_column_itd_thermodynamics, &
         config_calc_surface_temperature, &
         config_use_vertical_tracers

    real(kind=RKIND), pointer :: &
         config_dt

    call MPAS_pool_get_config(domain % configs, "config_use_column_package", config_use_column_package)

    if (config_use_column_package) then

       call MPAS_pool_get_config(domain % configs, "config_use_column_shortwave", config_use_column_shortwave)
       call MPAS_pool_get_config(domain % configs, "config_use_column_vertical_thermodynamics", &
                                                    config_use_column_vertical_thermodynamics)
       call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
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
       if (config_use_column_biogeochemistry) &
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
       call cice_column_update_state(domain, "thermodynamics", config_dt, config_dt)
       call mpas_timer_stop("Column predyn update state")

       !-----------------------------------------------------------------
       ! Separate vertical snow and ice tracers for advection
       !-----------------------------------------------------------------

       call mpas_timer_start("Column separate snow/ice tracers")
       if (config_use_vertical_tracers) &
            call column_separate_snow_ice_tracers(domain)
       call mpas_timer_stop("Column separate snow/ice tracers")

    endif

  end subroutine cice_column_predynamics_time_integration

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_dynamics_time_integration
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 6th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_column_dynamics_time_integration(domain, clock)

    type(domain_type), intent(inout) :: domain

    type(MPAS_clock_type), intent(in) :: clock

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_ridging, &
         config_use_vertical_tracers

    type(MPAS_pool_type), pointer :: &
         velocitySolver

    real(kind=RKIND), pointer :: &
         dynamicsTimeStep

    call MPAS_pool_get_config(domain % configs, "config_use_column_package", config_use_column_package)

    if (config_use_column_package) then

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
       call cice_column_update_state(domain, "transport", dynamicsTimeStep, 0.0_RKIND)
       call mpas_timer_stop("Column update state")

    endif

  end subroutine cice_column_dynamics_time_integration

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_postdynamics_time_integration
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 6th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_column_postdynamics_time_integration(domain, clock)

    type(domain_type), intent(inout) :: domain

    type(MPAS_clock_type), intent(in) :: clock

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_shortwave

    call MPAS_pool_get_config(domain % configs, "config_use_column_package", config_use_column_package)

    if (config_use_column_package) then

       call MPAS_pool_get_config(domain % configs, "config_use_column_shortwave", config_use_column_shortwave)

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
       call cice_column_coupling_prep(domain)
       call mpas_timer_stop("Column coupling prep")

    endif

  end subroutine cice_column_postdynamics_time_integration

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

    use ice_colpkg, only: &
         colpkg_step_therm1, &
         colpkg_clear_warnings, &
         colpkg_get_warnings

    use cice_constants, only: &
         cicePuny

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
         diagnostics

    ! configs
    real(kind=RKIND), pointer :: &
         config_dt

    logical, pointer :: &
         config_use_aerosols

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
         evaporativeWaterFlux, &
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
         freezeOnset, &
         oceanHeatFluxIceBottom

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
         sensibleHeatFluxCouple, &
         surfaceHeatFluxCouple, &
         surfaceConductiveFluxCouple, &
         snowThicknessChangeCategory

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
         iceBodyAerosol

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

    character(len=strKINDWarnings), dimension(:), allocatable :: &
         warnings

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

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)

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

       call MPAS_pool_get_array(ocean_fluxes, "oceanFreshWaterFlux", oceanFreshWaterFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanSaltFlux", oceanSaltFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanHeatFlux", oceanHeatFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanShortwaveFlux", oceanShortwaveFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanHeatFluxIceBottom", oceanHeatFluxIceBottom)

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

       ! loop over cells
       do iCell = 1, nCellsSolve

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

          abortFlag = .false.
          abortMessage = ""

          call colpkg_clear_warnings()
          call colpkg_step_therm1(&
               config_dt, &
               nCategories, &
               nIceLayers, &
               nSnowLayers, &
               nAerosols, &
               iceAreaCategoryInitial(:,iCell), &
               iceVolumeCategoryInitial(:,iCell), &
               snowVolumeCategoryInitial(:,iCell), &
               iceAreaCell(iCell), &
               iceAreaCategory(1,:,iCell), &
               iceVolumeCell(iCell), &
               iceVolumeCategory(1,:,iCell), &
               snowVolumeCell(iCell), &
               snowVolumeCategory(1,:,iCell), &
               uVelocity(iCell), &
               vVelocity(iCell), &
               surfaceTemperature(1,:,iCell), &
               snowEnthalpy(:,:,iCell), &
               iceEnthalpy(:,:,iCell), &
               iceSalinity(:,:,iCell), &
               levelIceArea(1,:,iCell), &
               levelIceVolume(1,:,iCell), &
               pondArea(1,:,iCell), &
               pondDepth(1,:,iCell), &
               pondLidThickness(1,:,iCell), &
               iceAge(1,:,iCell), &
               firstYearIceArea(1,:,iCell), &
               specificSnowAerosol(:,:,:), &
               specificIceAerosol(:,:,:), &
               uAirVelocity(iCell), &
               vAirVelocity(iCell), &
               windSpeed(iCell), &
               airLevelHeight(iCell), &
               airSpecificHumidity(iCell), &
               airDensity(iCell), &
               airTemperature(iCell), &
               atmosReferenceTemperature2m(iCell), &
               atmosReferenceHumidity2m(iCell), &
               atmosReferenceSpeed10m(iCell), &
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
               oceanHeatFluxIceBottom(iCell), &
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
               absorbedShortwaveSnowLayer(:,:,iCell), &
               absorbedShortwaveIceLayer(:,:,iCell), &
               longwaveDown(iCell), &
               solarZenithAngleCosine(iCell), &
               sensibleHeatFlux(iCell), &
               sensibleHeatFluxCategory(:,iCell), &
               latentHeatFlux(iCell), &
               latentHeatFluxCategory(:,iCell), &
               evaporativeWaterFlux(iCell), &
               oceanFreshWaterFlux(iCell), &
               oceanSaltFlux(iCell), &
               oceanHeatFlux(iCell), &
               oceanShortwaveFlux(iCell), &
               latentHeatFluxCouple(:,iCell), &
               sensibleHeatFluxCouple(:,iCell), &
               surfaceHeatFluxCouple(:,iCell), &
               surfaceConductiveFluxCouple(:,iCell), &
               atmosAerosolFlux(:,iCell), &
               oceanAerosolFlux(:,iCell), &
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
               snowThicknessChangeCategory(:,iCell), &
               frazilFormation(iCell), &
               northernHemisphereMask, &
               .not. northernHemisphereMask, &
               meltOnset(iCell), &
               freezeOnset(iCell), &
               dayOfYear, &
               abortFlag, &
               abortMessage)
          call colpkg_get_warnings(warnings)
          call column_write_warnings(warnings)

          ! code abort
          if (abortFlag) then
             COLUMN_ERROR_WRITE("column_vertical_thermodynamics")
             COLUMN_ERROR_WRITE(trim(abortMessage))
             write(abortLocation,*) indexToCellID(iCell)
             call MPAS_dmpar_global_abort("MPAS-seaice: iCell "//trim(abortLocation)//": "//trim(abortMessage))
          endif

          ! aerosol
          if (config_use_aerosols) then

             do iCategory = 1, nCategories
                do iAerosol = 1, nAerosols

                   if (snowVolumeCategory(1,iCategory,iCell) > cicePuny) &
                        specificSnowAerosol(iAerosol, :, iCategory) = &
                        specificSnowAerosol(iAerosol, :, iCategory) / snowVolumeCategory(1,iCategory,iCell)

                   if (iceVolumeCategory(1,iCategory,iCell) > cicePuny) &
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
!> \author Adrian K. Turner, LANL
!> \date 21th January 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine column_itd_thermodynamics(domain, clock)

    use ice_colpkg, only: &
         colpkg_step_therm2, &
         colpkg_get_warnings, &
         colpkg_clear_warnings

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
         config_use_column_biogeochemistry

    ! dimensions
    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nIceLayers, &
         nSnowLayers, &
         nAerosols, &
         nBioLayers

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
         freezeOnset, &
         categoryThicknessLimits, &
         biologyGrid, &
         verticalGrid, &
         interfaceBiologyGrid, &
         zSalinityFlux, &
         frazilGrowthDiagnostic

    real(kind=RKIND), dimension(:,:), pointer :: &
         iceAreaCategoryInitial, &
         iceVolumeCategoryInitial, &
         oceanAerosolFlux, &
         oceanBioFluxes, &
         oceanBioConcentrations, &
         initialSalinityProfile

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    integer, dimension(:,:), pointer :: &
         newlyFormedIce

    integer, dimension(:), pointer :: &
         indexToCellID

    ! local
    integer :: &
         iCell, &
         iCategory, &
         iBioTracers, &
         iBioData

    logical, dimension(:), allocatable :: &
         newlyFormedIceLogical

    real(kind=RKIND), dimension(:), allocatable :: &
         oceanBioConcentrationsUsed

    logical :: &
         abortFlag, &
         setGetPhysicsTracers, &
         setGetBGCTracers

    character(len=strKIND) :: &
         abortMessage, &
         abortLocation

    real(kind=RKIND) :: &
         dayOfYear

    character(len=strKINDWarnings), dimension(:), allocatable :: &
         warnings

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

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_update_ocean_fluxes", config_update_ocean_fluxes)
       call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)
       call MPAS_pool_get_dimension(mesh, "nAerosols", nAerosols)
       call MPAS_pool_get_dimension(block % dimensions, "nBioLayers", nBioLayers)

       call MPAS_pool_get_array(mesh, "indexToCellID", indexToCellID)

       call MPAS_pool_get_array(icestate, "iceAreaCategoryInitial", iceAreaCategoryInitial)
       call MPAS_pool_get_array(icestate, "iceVolumeCategoryInitial", iceVolumeCategoryInitial)
       call MPAS_pool_get_array(icestate, "openWaterArea", openWaterArea)

       call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

       call MPAS_pool_get_array(atmos_coupling, "rainfallRate", rainfallRate)

       call MPAS_pool_get_array(ocean_coupling, "freezingMeltingPotential", freezingMeltingPotential)
       call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)
       call MPAS_pool_get_array(ocean_coupling, "seaSurfaceSalinity", seaSurfaceSalinity)

       call MPAS_pool_get_array(ocean_fluxes, "oceanFreshWaterFlux", oceanFreshWaterFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanSaltFlux", oceanSaltFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanHeatFlux", oceanHeatFlux)

       call MPAS_pool_get_array(melt_growth_rates, "lateralIceMeltFraction", lateralIceMeltFraction)
       call MPAS_pool_get_array(melt_growth_rates, "lateralIceMelt", lateralIceMelt)
       call MPAS_pool_get_array(melt_growth_rates, "frazilFormation", frazilFormation)
       call MPAS_pool_get_array(melt_growth_rates, "frazilGrowthDiagnostic", frazilGrowthDiagnostic)

       call MPAS_pool_get_array(ponds, "pondFreshWaterFlux", pondFreshWaterFlux)

       call MPAS_pool_get_array(aerosols, "oceanAerosolFlux", oceanAerosolFlux)

       call MPAS_pool_get_array(biogeochemistry, "newlyFormedIce", newlyFormedIce)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioFluxes", oceanBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioConcentrations", oceanBioConcentrations)
       call MPAS_pool_get_array(biogeochemistry, "biologyGrid", biologyGrid)
       call MPAS_pool_get_array(biogeochemistry, "verticalGrid", verticalGrid)
       call MPAS_pool_get_array(biogeochemistry, "interfaceBiologyGrid", interfaceBiologyGrid)
       call MPAS_pool_get_array(biogeochemistry, "zSalinityFlux", zSalinityFlux)

       call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)
       call MPAS_pool_get_array(initial, "categoryThicknessLimits", categoryThicknessLimits)

       call MPAS_pool_get_array(diagnostics, "freezeOnset", freezeOnset)

       ! newly formed ice
       allocate(newlyFormedIceLogical(nCategories))
       allocate(oceanBioConcentrationsUsed(ciceTracerObject % nBioTracers))

       setGetPhysicsTracers = .true.
       setGetBGCTracers     = config_use_column_biogeochemistry

       do iCell = 1, nCellsSolve

          ! newly formed ice
          do iCategory = 1, nCategories
             newlyFormedIceLogical(iCategory) = (newlyFormedIce(iCategory,iCell) == 1)
          enddo ! iCategory

          ! read the required ocean concentration fields into the allocated array
          do iBioTracers = 1, ciceTracerObject % nBioTracers
             iBioData = ciceTracerObject % index_LayerIndexToDataArray(iBioTracers)
             oceanBioConcentrationsUsed(iBioTracers) = oceanBioConcentrations(iBioData, iCell)
          enddo ! iBioTracers

          ! set the category tracer array
          call set_cice_tracer_array_category(block, ciceTracerObject, iCell, setGetPhysicsTracers, setGetBGCTracers)

          abortFlag = .false.
          abortMessage = ""

          call colpkg_clear_warnings()
          call colpkg_step_therm2(&
               config_dt, &
               nCategories, &
               nAerosols, &
               ciceTracerObject % nBioTracersLayer, &
               nIcelayers, &
               nSnowLayers, &
               categoryThicknessLimits(:), &
               nBioLayers, &
               iceAreaCategory(1,:,iCell), &
               iceVolumeCategory(1,:,iCell), &
               snowVolumeCategory(1,:,iCell), &
               iceAreaCategoryInitial(:,iCell), &
               iceVolumeCategoryInitial(:,iCell), &
               ciceTracerObject % tracerArrayCategory, & !trcrn
               openWaterArea(iCell), &
               iceAreaCell(iCell), &
               ciceTracerObject % parentIndex, & !trcr_depend
               ciceTracerObject % firstAncestorMask, & !trcr_base
               ciceTracerObject % ancestorNumber, & !n_trcr_strata
               ciceTracerObject % ancestorIndices, & !nt_strata
               seaFreezingTemperature(iCell), &
               seaSurfaceSalinity(iCell), &
               initialSalinityProfile(:,iCell), &
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
               biologyGrid(:), &
               verticalGrid(:), &
               interfaceBiologyGrid(:), &
               oceanAerosolFlux(:,iCell), &
               newlyFormedIceLogical(:), &
               zSalinityFlux(iCell), &
               oceanBioFluxes(:,iCell), &
               oceanBioConcentrationsUsed(:), &
               abortFlag, &
               abortMessage, &
               frazilGrowthDiagnostic(iCell), &
               freezeOnset(iCell), &
               dayOfYear)
          call colpkg_get_warnings(warnings)
          call column_write_warnings(warnings)

          ! update
          do iCategory = 1, nCategories
             newlyFormedIce(iCategory,iCell) = 0
             if (newlyFormedIceLogical(iCategory)) newlyFormedIce(iCategory,iCell) = 1
          enddo ! iCategory

          ! get category tracer array
          call get_cice_tracer_array_category(block, ciceTracerObject, iCell, setGetPhysicsTracers, setGetBGCTracers)

          ! code abort
          if (abortFlag) then
             COLUMN_ERROR_WRITE("column_itd_thermodynamics")
             COLUMN_ERROR_WRITE(trim(abortMessage))
             write(abortLocation,*) indexToCellID(iCell)
             call MPAS_dmpar_global_abort("MPAS-seaice: iCell "//trim(abortLocation)//": "//trim(abortMessage))
          endif

       enddo ! iCell

       ! newly formed ice
       deallocate(newlyFormedIceLogical)
       deallocate(oceanBioConcentrationsUsed)

       block => block % next
    end do

  end subroutine column_itd_thermodynamics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_prep_radiation
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 21th January 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine column_prep_radiation(domain)

    use ice_colpkg, only: colpkg_prep_radiation

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
         penetratingShortwaveFlux

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
       call MPAS_pool_get_array(shortwave, "shortwaveLayerPenetration", shortwaveLayerPenetration)
       call MPAS_pool_get_array(shortwave, "absorbedShortwaveSnowLayer", absorbedShortwaveSnowLayer)
       call MPAS_pool_get_array(shortwave, "absorbedShortwaveIceLayer", absorbedShortwaveIceLayer)

       do iCell = 1, nCellsSolve

          call colpkg_prep_radiation(&
               nCategories, &
               nIceLayers, &
               nSnowLayers, &
               iceAreaCell(iCell), &
               iceAreaCategory(1,:,iCell), &
               shortwaveVisibleDirectDown(iCell), &
               shortwaveVisibleDiffuseDown(iCell), &
               shortwaveIRDirectDown(iCell), &
               shortwaveIRDiffuseDown(iCell), &
               albedoVisibleDirectArea(iCell), &
               albedoVisibleDiffuseArea(iCell), &
               albedoIRDirectArea(iCell), &
               albedoIRDiffuseArea(iCell), &
               shortwaveScalingFactor(iCell), &
               surfaceShortwaveFlux(:,iCell), &
               interiorShortwaveFlux(:,iCell), &
               penetratingShortwaveFlux(:,iCell), &
               shortwaveLayerPenetration(:,:,iCell), &
               absorbedShortwaveSnowLayer(:,:,iCell), &
               absorbedShortwaveIceLayer(:,:,iCell))

       enddo ! iCell

       block => block % next
    end do

  end subroutine column_prep_radiation

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_radiation
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 21th January 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine column_radiation(domain, clock, lInitialization)

    use ice_colpkg, only: &
         colpkg_step_radiation, &
         colpkg_get_warnings, &
         colpkg_clear_warnings

    use cice_constants, only: pii

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
         aerosols, &
         biogeochemistry

    ! configs
    real(kind=RKIND), pointer :: &
         config_dt

    logical, pointer :: &
         config_use_shortwave_bioabsorption, &
         config_use_brine, &
         config_use_modal_aerosols, &
         config_use_zaerosols, &
         config_use_column_biogeochemistry

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
         snowfallRate, &
         verticalShortwaveGrid, &
         verticalGrid

    real(kind=RKIND), dimension(:,:), pointer :: &
         surfaceShortwaveFlux, &
         interiorShortwaveFlux, &
         penetratingShortwaveFlux, &
         bareIceAlbedoCategory, &
         snowAlbedoCategory, &
         pondAlbedoCategory, &
         effectivePondAreaCategory, &
         pondSnowDepthDifference, &
         pondLidMeltFluxFraction, &
         aerosolMassExtinctionCrossSection, &
         aerosolSingleScatterAlbedo, &
         aerosolAsymmetryParameter, &
         modalMassExtinctionCrossSection, &
         modalSingleScatterAlbedo, &
         modalAsymmetryParameter, &
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
         modalBCabsorptionParameter, &
         bioTracerShortwave

    real(kind=RKIND), pointer :: &
         dayOfNextShortwaveCalculation ! needed for CESM like coupled simulations

    character(len=strKIND), pointer :: &
         config_calendar_type

    character(len=strKIND) :: &
         calendarType ! needed for CESM like coupled simulations

    ! aerosols array
    real(kind=RKIND), dimension(:,:), allocatable :: &
         aerosolsArray

    ! local
    integer :: &
         iCell, &
         iCategory, &
         iAerosol, &
         iTracer

    integer, dimension(:), allocatable :: &
         index_shortwaveAerosol

    real(kind=RKIND) :: &
         dayOfYear, &
         lonCellColumn

    integer :: &
         secondsIntoDay, &
         daysInYear

    logical :: &
         setGetPhysicsTracers, &
         setGetBGCTracers

    character(len=strKINDWarnings), dimension(:), allocatable :: &
         warnings

    ! day of year
    call get_day_of_year(clock, dayOfYear)

    ! seconds into day
    call get_seconds_into_day(clock, secondsIntoDay)

    ! get days in year
    call get_days_in_year(domain, clock, daysInYear)

    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)
    call MPAS_pool_get_config(domain % configs, "config_use_shortwave_bioabsorption", config_use_shortwave_bioabsorption)
    call MPAS_pool_get_config(domain % configs, "config_use_modal_aerosols",config_use_modal_aerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_zaerosols",config_use_zaerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry",config_use_column_biogeochemistry)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "ponds", ponds)
       call MPAS_pool_get_subpool(block % structs, "aerosols", aerosols)
       call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)
       call MPAS_pool_get_dimension(mesh, "nAerosols", nAerosols)
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

       call MPAS_pool_get_array(aerosols, "aerosolMassExtinctionCrossSection", aerosolMassExtinctionCrossSection)
       call MPAS_pool_get_array(aerosols, "aerosolSingleScatterAlbedo", aerosolSingleScatterAlbedo)
       call MPAS_pool_get_array(aerosols, "aerosolAsymmetryParameter", aerosolAsymmetryParameter)
       call MPAS_pool_get_array(aerosols, "modalMassExtinctionCrossSection", modalMassExtinctionCrossSection)
       call MPAS_pool_get_array(aerosols, "modalSingleScatterAlbedo", modalSingleScatterAlbedo)
       call MPAS_pool_get_array(aerosols, "modalAsymmetryParameter", modalAsymmetryParameter)
       call MPAS_pool_get_array(aerosols, "modalBCabsorptionParameter", modalBCabsorptionParameter)

       call MPAS_pool_get_array(biogeochemistry, "bioTracerShortwave", bioTracerShortwave)
       call MPAS_pool_get_array(biogeochemistry, "verticalShortwaveGrid", verticalShortwaveGrid)
       call MPAS_pool_get_array(biogeochemistry, "verticalGrid", verticalGrid)

       ! calendar type
       call MPAS_pool_get_config(block % configs, "config_calendar_type", config_calendar_type)
       if (trim(config_calendar_type) == "gregorian") then
          calendarType = "GREGORIAN"
       else
          calendarType = "GREGORIAN_NOLEAP"
       endif

       ! aerosols array
       allocate(aerosolsArray(4*nAerosols,nCategories))
       allocate(index_shortwaveAerosol(maxAerosolType))
       if (.not. config_use_column_biogeochemistry) then
          index_shortwaveAerosol(1:maxAerosolType) = 1
       else
          do iAerosol = 1, maxAerosolType
               index_shortwaveAerosol(iAerosol) =  ciceTracerObject % index_verticalAerosolsConcShortwave(iAerosol)
          enddo
       endif

       setGetPhysicsTracers = .true.
       setGetBGCTracers     = config_use_column_biogeochemistry

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

          lonCellColumn = lonCell(iCell)
          if (lonCellColumn > pii) lonCellColumn = lonCellColumn - 2.0_RKIND * pii

          ! set the category tracer array
          call set_cice_tracer_array_category(block, ciceTracerObject, iCell, setGetPhysicsTracers, setGetBGCTracers)

          call colpkg_clear_warnings()
          call colpkg_step_radiation(&
               config_dt, &
               nCategories, &
               nAlgae, &
               config_use_zaerosols, &
               nBioLayers, &
               ciceTracerObject % nTracers, &
               ciceTracerObject % nBioTracers, &
               ciceTracerObject % nBioTracersShortwave, &
               nIceLayers, &
               nSnowLayers, &
               nAerosols, &
               nzAerosols, &
               config_use_shortwave_bioabsorption, &
               ciceTracerObject % index_chlorophyllShortwave, &
               index_shortwaveAerosol, &
               verticalShortwaveGrid(:), &
               verticalGrid(:), &
               brineFraction(1,:,iCell), &
               iceAreaCategory(1,:,iCell), &
               iceVolumeCategory(1,:,iCell), &
               snowVolumeCategory(1,:,iCell), &
               surfaceTemperature(1,:,iCell), &
               levelIceArea(1,:,iCell), &
               pondArea(1,:,iCell), &
               pondDepth(1,:,iCell), &
               pondLidThickness(1,:,iCell), &
               aerosolsArray, &
               bioTracerShortwave(:,:,iCell), &
               ciceTracerObject % tracerArrayCategory(:,:), &
               latCell(iCell), &
               lonCellColumn, &
               calendarType, &
               daysInYear, &
               dayOfNextShortwaveCalculation, &
               dayOfYear, &
               secondsIntoDay, &
               aerosolMassExtinctionCrossSection(:,:), &
               aerosolSingleScatterAlbedo(:,:), &
               aerosolAsymmetryParameter(:,:), &
               modalMassExtinctionCrossSection(:,:), &
               modalSingleScatterAlbedo(:,:), &
               modalAsymmetryParameter(:,:), &
               modalBCabsorptionParameter(:,:,:), &
               config_use_modal_aerosols, &
               shortwaveVisibleDirectDown(iCell), &
               shortwaveVisibleDiffuseDown(iCell), &
               shortwaveIRDirectDown(iCell), &
               shortwaveIRDiffuseDown(iCell), &
               solarZenithAngleCosine(iCell), &
               snowfallRate(iCell), &
               albedoVisibleDirectCategory(:,iCell), &
               albedoVisibleDiffuseCategory(:,iCell), &
               albedoIRDirectCategory(:,iCell), &
               albedoIRDiffuseCategory(:,iCell), &
               surfaceShortwaveFlux(:,iCell), &
               interiorShortwaveFlux(:,iCell), &
               penetratingShortwaveFlux(:,iCell), &
               shortwaveLayerPenetration(:,:,iCell), &
               absorbedShortwaveSnowLayer(:,:,iCell), &
               absorbedShortwaveIceLayer(:,:,iCell), &
               bareIceAlbedoCategory(:,iCell), &
               snowAlbedoCategory(:,iCell), &
               pondAlbedoCategory(:,iCell), &
               effectivePondAreaCategory(:,iCell), &
               snowFractionCategory(:,iCell), &
               pondSnowDepthDifference(:,iCell), &
               pondLidMeltFluxFraction(:,iCell), &
               .false., &
               lInitialization)
          call colpkg_get_warnings(warnings)
          call column_write_warnings(warnings)

          ! set the category tracer array
          call get_cice_tracer_array_category(block, ciceTracerObject, iCell, setGetPhysicsTracers, setGetBGCTracers)

       enddo ! iCell

       ! aerosols array
       deallocate(aerosolsArray)
       deallocate(index_shortwaveAerosol)

       block => block % next
    end do

  end subroutine column_radiation

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_ridging
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 21th January 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine column_ridging(domain)

    use ice_colpkg, only: &
         colpkg_step_ridge, &
         colpkg_get_warnings, &
         colpkg_clear_warnings

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         icestate, &
         tracers, &
         tracers_aggregate, &
         ponds, &
         ocean_fluxes, &
         ridging, &
         aerosols, &
         biogeochemistry, &
         initial, &
         velocity_solver

    ! configs
    logical, pointer :: &
         config_use_column_biogeochemistry

    real(kind=RKIND), pointer :: &
         config_dt

    integer, pointer :: &
         config_dynamics_subcycle_number

    ! dimensions
    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nIceLayers, &
         nSnowLayers, &
         nAerosols, &
         nBioLayers

    ! variables
    real(kind=RKIND), dimension(:), pointer :: &
         pondFreshWaterFlux, &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         iceAreaCell, &
         ridgeConvergence, &
         ridgeShear, &
         openWaterArea, &
         areaLossRidge, &
         areaGainRidge, &
         iceVolumeRidged, &
         openingRateRidge, &
         categoryThicknessLimits, &
         zSalinityFlux

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
         iCategory

    logical, dimension(:), allocatable :: &
         newlyFormedIceLogical

    logical :: &
         abortFlag, &
         setGetPhysicsTracers, &
         setGetBGCTracers

    character(len=strKIND) :: &
         abortMessage, &
         abortLocation

    character(len=strKINDWarnings), dimension(:), allocatable :: &
         warnings

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

       call MPAS_pool_get_config(block % configs, "config_dynamics_subcycle_number", config_dynamics_subcycle_number)
       call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

       call MPAS_pool_get_array(velocity_solver, "dynamicsTimeStep", dynamicsTimeStep)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)
       call MPAS_pool_get_dimension(mesh, "nAerosols", nAerosols)
       call MPAS_pool_get_dimension(block % dimensions, "nBioLayers", nBioLayers)

       call MPAS_pool_get_array(mesh, "indexToCellID", indexToCellID)

       call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)

       call MPAS_pool_get_array(icestate, "openWaterArea", openWaterArea)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

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
       call MPAS_pool_get_array(biogeochemistry, "zSalinityFlux", zSalinityFlux)

       call MPAS_pool_get_array(initial, "categoryThicknessLimits", categoryThicknessLimits)

       ! newly formed ice
       allocate(newlyFormedIceLogical(nCategories))

       setGetPhysicsTracers = .true.
       setGetBGCTracers     = config_use_column_biogeochemistry

       do iCell = 1, nCellsSolve

          ! newly formed ice
          do iCategory = 1, nCategories
             newlyFormedIceLogical(iCategory) = (newlyFormedIce(iCategory,iCell) == 1)
          enddo ! iCategory

          ! set the category tracer array
          call set_cice_tracer_array_category(block, ciceTracerObject, iCell, setGetPhysicsTracers, setGetBGCTracers)

          abortFlag = .false.
          abortMessage = ""

          call colpkg_clear_warnings()
          call colpkg_step_ridge(&
               dynamicsTimeStep, &
               config_dynamics_subcycle_number, &
               nIceLayers, &
               nSnowLayers, &
               nBioLayers, &
               nCategories, &
               categoryThicknessLimits, &
               ridgeConvergence(iCell), &
               ridgeShear(iCell), &
               iceAreaCategory(1,:,iCell), &
               ciceTracerObject % tracerArrayCategory, & ! trcrn
               iceVolumeCategory(1,:,iCell), &
               snowVolumeCategory(1,:,iCell), &
               openWaterArea(iCell), &
               ciceTracerObject % parentIndex, & ! trcr_depend
               ciceTracerObject % firstAncestorMask, & ! trcr_base
               ciceTracerObject % ancestorNumber, & ! n_trcr_strata
               ciceTracerObject % ancestorIndices, & ! nt_strata
               areaLossRidge(iCell), &
               areaGainRidge(iCell), &
               iceVolumeRidged(iCell), &
               openingRateRidge(iCell), &
               pondFreshWaterFlux(iCell), &
               oceanFreshWaterFlux(iCell), &
               oceanHeatFlux(iCell), &
               nAerosols, &
               oceanAerosolFlux(:,iCell), &
               ridgeParticipationFunction(:,iCell), &
               ratioRidgeThicknessToIce(:,iCell), &
               fractionNewRidgeArea(:,iCell), &
               fractionNewRidgeVolume(:,iCell), &
               areaLossRidgeCategory(:,iCell), &
               areaGainRidgeCategory(:,iCell), &
               iceVolumeRidgedCategory(:,iCell), &
               raftingIceArea(:,iCell), &
               raftingIceVolume(:,iCell), &
               iceAreaCell(iCell), &
               oceanSaltFlux(iCell), &
               newlyFormedIceLogical(:), &
               zSalinityFlux(iCell), &
               oceanBioFluxes(:,iCell), &
               abortFlag, &
               abortMessage)
          call colpkg_get_warnings(warnings)
          call column_write_warnings(warnings)

          ! update
          do iCategory = 1, nCategories
             newlyFormedIce(iCategory,iCell) = 0
             if (newlyFormedIceLogical(iCategory)) newlyFormedIce(iCategory,iCell) = 1
          enddo ! iCategory

          ! get category tracer array
          call get_cice_tracer_array_category(block, ciceTracerObject, iCell, setGetPhysicsTracers, setGetBGCTracers)

          ! code abort
          if (abortFlag) then
             COLUMN_ERROR_WRITE("column_ridging")
             COLUMN_ERROR_WRITE(trim(abortMessage))
             write(abortLocation,*) indexToCellID(iCell)
             call MPAS_dmpar_global_abort("MPAS-seaice: iCell "//trim(abortLocation)//": "//trim(abortMessage))
          endif

       enddo ! iCell

       ! newly formed ice
       deallocate(newlyFormedIceLogical)

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

    use ice_colpkg, only: &
         colpkg_biogeochemistry, &
         colpkg_init_OceanConcArray, &
         colpkg_get_warnings, &
         colpkg_clear_warnings

    use cice_constants, only: &
         cicePuny

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         biogeochemistry, &
         icestate, &
         tracers, &
         shortwave, &
         melt_growth_rates, &
         ocean_coupling, &
         atmos_coupling, &
         initial

    ! configs
    real(kind=RKIND), pointer :: &
         config_dt

    logical, pointer :: &
         config_use_brine, &
         config_use_skeletal_biochemistry, &
         config_use_column_biogeochemistry, &
         config_use_zaerosols

    ! dimensions
    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nIceLayers, &
         nSnowLayers, &
         nzAerosols, &
         nBioLayers, &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON, &
         nParticulateIron, &
         nDissolvedIron, &
         nZBGCTracers, &
         maxAlgaeType, &
         maxDOCType, &
         maxDICType, &
         maxDONType, &
         maxIronType, &
         maxBCType, &
         maxDustType, &
         maxAerosolType

    ! variables

    real(kind=RKIND), dimension(:), pointer :: &
         rayleighCriteriaReal, &
         netNitrateUptake, &
         netAmmoniumUptake, &
         totalVerticalSalinity, &
         netSpecificAlgalGrowthRate, &
         primaryProduction, &
         netBrineHeight, &
         biologyGrid, &
         interfaceBiologyGrid, &
         interfaceGrid, &
         verticalGrid, &
         seaSurfaceTemperature, &
         seaSurfaceSalinity, &
         snowfallRate, &
         zSalinityFlux, &
         zSalinityGDFlux, &
         oceanMixedLayerDepth, &
         totalSkeletalAlgae, &
         oceanNitrateConc, &
         oceanSilicateConc, &
         oceanAmmoniumConc, &
         oceanDMSConc, &
         oceanDMSPConc, &
         oceanHumicsConc

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
         totalVerticalBiologyIce, &
         totalVerticalBiologySnow, &
         penetratingShortwaveFlux, &
         zSalinityIceDensity, &
         basalIceMeltCategory, &
         surfaceIceMeltCategory, &
         congelationCategory, &
         snowiceFormationCategory, &
         snowMeltCategory, &
         initialSalinityProfile, &
         atmosBioFluxes, &
         atmosBlackCarbonFlux, &
         atmosDustFlux, &
         oceanBioFluxes, &
         oceanAlgaeConc, &
         oceanDOCConc, &
         oceanDICConc, &
         oceanDONConc, &
         oceanParticulateIronConc, &
         oceanDissolvedIronConc, &
         oceanZAerosolConc

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
         skeletalAlgaeConc

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
         indexj

    real(kind=RKIND), dimension(:), allocatable :: &
         oceanBioConcentrationsUsed

    logical, dimension(:), allocatable :: &
         newlyFormedIceLogical

    logical :: &
         abortFlag, &
         rayleighCriteria, &
         setGetPhysicsTracers, &
         setGetBGCTracers

    character(len=strKIND) :: &
         abortMessage, &
         abortLocation

    character(len=strKINDWarnings), dimension(:), allocatable :: &
         warnings

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
       call MPAS_pool_get_subpool(block % structs, "melt_growth_rates", melt_growth_rates)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
       call MPAS_pool_get_subpool(block % structs, "initial", initial)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)
       call MPAS_pool_get_dimension(mesh, "nzAerosols", nzAerosols)
       call MPAS_pool_get_dimension(mesh, "nBioLayers", nBioLayers)
       call MPAS_pool_get_dimension(mesh, "nAlgae", nAlgae)
       call MPAS_pool_get_dimension(mesh, "nDOC", nDOC)
       call MPAS_pool_get_dimension(mesh, "nDIC", nDIC)
       call MPAS_pool_get_dimension(mesh, "nDON", nDON)
       call MPAS_pool_get_dimension(mesh, "nParticulateIron", nParticulateIron)
       call MPAS_pool_get_dimension(mesh, "nDissolvedIron", nDissolvedIron)
       call MPAS_pool_get_dimension(mesh, "nZBGCTracers", nZBGCTracers)
       call MPAS_pool_get_dimension(mesh, "maxAlgaeType", maxAlgaeType)
       call MPAS_pool_get_dimension(mesh, "maxDOCType", maxDOCType)
       call MPAS_pool_get_dimension(mesh, "maxDICType", maxDICType)
       call MPAS_pool_get_dimension(mesh, "maxDONType", maxDONType)
       call MPAS_pool_get_dimension(mesh, "maxAerosolType", maxAerosolType)
       call MPAS_pool_get_dimension(mesh, "maxIronType", maxIronType)
       call MPAS_pool_get_dimension(mesh, "maxBCType", maxBCType)
       call MPAS_pool_get_dimension(mesh, "maxDustType", maxDustType)

       call MPAS_pool_get_dimension(mesh, "indexToCellID", indexToCellID)

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_use_brine", config_use_brine)
       call MPAS_pool_get_config(block % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
       call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
       call MPAS_pool_get_config(block % configs, "config_use_zaerosols",config_use_zaerosols)

       call MPAS_pool_get_array(biogeochemistry, "newlyFormedIce", newlyFormedIce)
       call MPAS_pool_get_array(biogeochemistry, "netNitrateUptake", netNitrateUptake)
       call MPAS_pool_get_array(biogeochemistry, "netAmmoniumUptake", netAmmoniumUptake)
       call MPAS_pool_get_array(biogeochemistry, "totalVerticalSalinity", totalVerticalSalinity)
       call MPAS_pool_get_array(biogeochemistry, "netSpecificAlgalGrowthRate", netSpecificAlgalGrowthRate)
       call MPAS_pool_get_array(biogeochemistry, "primaryProduction", primaryProduction)
       call MPAS_pool_get_array(biogeochemistry, "netBrineHeight", netBrineHeight)
       call MPAS_pool_get_array(biogeochemistry, "brineBottomChange", brineBottomChange)
       call MPAS_pool_get_array(biogeochemistry, "brineTopChange", brineTopChange)
       call MPAS_pool_get_array(biogeochemistry, "bioPorosity", bioPorosity)
       call MPAS_pool_get_array(biogeochemistry, "rayleighCriteriaReal", rayleighCriteriaReal)
       call MPAS_pool_get_array(biogeochemistry, "biologyGrid", biologyGrid)
       call MPAS_pool_get_array(biogeochemistry, "interfaceBiologyGrid", interfaceBiologyGrid)
       call MPAS_pool_get_array(biogeochemistry, "interfaceGrid", interfaceGrid)
       call MPAS_pool_get_array(biogeochemistry, "verticalGrid", verticalGrid)
       call MPAS_pool_get_array(biogeochemistry, "bioDiffusivity", bioDiffusivity)
       call MPAS_pool_get_array(biogeochemistry, "bioPermeability", bioPermeability)
       call MPAS_pool_get_array(biogeochemistry, "bioShortwaveFlux", bioShortwaveFlux)
       call MPAS_pool_get_array(biogeochemistry, "darcyVelocityBio", darcyVelocityBio)
       call MPAS_pool_get_array(biogeochemistry, "snowIceBioFluxes", snowIceBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "atmosIceBioFluxes", atmosIceBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioConcentrations", oceanBioConcentrations)
       call MPAS_pool_get_array(biogeochemistry, "totalVerticalBiologyIce", totalVerticalBiologyIce)
       call MPAS_pool_get_array(biogeochemistry, "totalVerticalBiologySnow", totalVerticalBiologySnow)
       call MPAS_pool_get_array(biogeochemistry, "zSalinityIceDensity", zSalinityIceDensity)
       call MPAS_pool_get_array(biogeochemistry, "zSalinityFlux", zSalinityFlux)
       call MPAS_pool_get_array(biogeochemistry, "zSalinityGDFlux", zSalinityGDFlux)
       call MPAS_pool_get_array(biogeochemistry, "atmosBioFluxes", atmosBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "atmosBlackCarbonFlux", atmosBlackCarbonFlux)
       call MPAS_pool_get_array(biogeochemistry, "atmosDustFlux", atmosDustFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioFluxes", oceanBioFluxes)
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

       call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTemperature", seaSurfaceTemperature)
       call MPAS_pool_get_array(ocean_coupling, "seaSurfaceSalinity", seaSurfaceSalinity)
       call MPAS_pool_get_array(ocean_coupling, "oceanMixedLayerDepth", oceanMixedLayerDepth)

       call MPAS_pool_get_array(atmos_coupling, "snowfallRate", snowfallRate)

       call MPAS_pool_get_array(icestate, "iceAreaCategoryInitial", iceAreaCategoryInitial)
       call MPAS_pool_get_array(icestate, "iceVolumeCategoryInitial", iceVolumeCategoryInitial)
       call MPAS_pool_get_array(icestate, "snowVolumeCategoryInitial", snowVolumeCategoryInitial)
       call MPAS_pool_get_array(icestate, "iceThicknessCategoryInitial", iceThicknessCategoryInitial)

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

       ! newly formed ice
       allocate(newlyFormedIceLogical(nCategories))
       allocate(oceanBioConcentrationsUsed(ciceTracerObject % nBioTracers))

       setGetPhysicsTracers = .true.
       setGetBGCTracers     = config_use_column_biogeochemistry

       do iCell = 1, nCellsSolve

          ! newly formed ice
          do iCategory = 1, nCategories
             newlyFormedIceLogical(iCategory) = (newlyFormedIce(iCategory,iCell) == 1)
          enddo ! iCategory
          rayleighCriteria = (rayleighCriteriaReal(iCell) > 0.5_RKIND)

          !update ocean concentrations fields and atmospheric fluxes into allocated array
#ifdef coupled
          call colpkg_init_OceanConcArray(&
               nZBGCTracers, &
               maxAlgaeType, &
               maxDONType, &
               maxDOCType, &
               maxDICType, &
               maxAerosolType, &
               maxIronType, &
               oceanNitrateConc(iCell), &
               oceanAmmoniumConc(iCell), &
               oceanSilicateConc(iCell),&
               oceanDMSPConc(iCell), &
               oceanDMSConc(iCell), &
               oceanAlgaeConc(:,iCell), &
               oceanDOCConc(:,iCell), &
               oceanDONConc(:,iCell), &
               oceanDICConc(:,iCell), &
               oceanDissolvedIronConc(:,iCell), &
               oceanParticulateIronConc(:,iCell), &
               oceanZAerosolConc(:,iCell), &
               oceanBioConcentrations(:,iCell), &
               oceanHumicsConc(iCell))
#else
          do iBioTracers = 1, maxBCType
             atmosBlackCarbonFlux(iBioTracers,iCell) =  1.e-12_RKIND
          enddo
          do iBioTracers = 1, maxDustType
             atmosDustFlux(iBioTracers,iCell) =  1.e-13_RKIND
          enddo
#endif

          if (config_use_zaerosols) then
             indexj = ciceTracerObject % index_verticalAerosolsConcLayer(1)
             do iBioTracers = 1, maxBCType
                atmosBioFluxes(indexj -1 + iBioTracers,iCell) = atmosBlackCarbonFlux(iBioTracers,iCell)
             enddo
             do iBioTracers = maxBCType + 1, nzAerosols
                atmosBioFluxes(indexj -1 + iBioTracers, iCell) = atmosDustFlux(iBioTracers-maxBCType,iCell)
             enddo
          endif

          do iBioTracers = 1, ciceTracerObject % nBioTracers
             iBioData = ciceTracerObject % index_LayerIndexToDataArray(iBioTracers)
             oceanBioConcentrationsUsed(iBioTracers) = oceanBioConcentrations(iBioData,iCell)
          enddo ! iBioTracers

          abortFlag = .false.

          call set_cice_tracer_array_category(block, ciceTracerObject, iCell, setGetPhysicsTracers, setGetBGCTracers)

          call colpkg_clear_warnings()
          call colpkg_biogeochemistry(&
               config_dt, &
               ciceTracerObject % nTracers, &
               ciceTracerObject % nBioTracers, &
               netNitrateUptake(iCell), &
               netAmmoniumUptake(iCell), &
               bioDiffusivity(:,:,iCell), &
               bioPermeability(:,:,iCell), &
               bioShortwaveFlux(:,:,iCell), &
               totalVerticalSalinity(iCell), &
               darcyVelocityBio(:,iCell), &
               netSpecificAlgalGrowthRate(iCell), &
               primaryProduction(iCell), &
               netBrineHeight(iCell), &
               brineBottomChange(:,iCell), &
               brineTopChange(:,iCell), &
               verticalNitrogenLosses(:,:,iCell), &
               snowIceBioFluxes(:,iCell), &
               atmosIceBioFluxes(:,iCell), &
               oceanBioConcentrationsUsed(:), &
               newlyFormedIceLogical(:), &
               shortwaveLayerPenetration(:,:,iCell), &
               bioPorosity(:,:,iCell), &
               bioTemperature(:,:,iCell), &
               totalVerticalBiologyIce(:,iCell), &
               totalVerticalBiologySnow(:,iCell), &
               penetratingShortwaveFlux(:,iCell), &
               rayleighCriteria, &
               zSalinityIceDensity(:,iCell), &
               zSalinityFlux(iCell), &
               zSalinityGDFlux(iCell), &
               biologyGrid, &
               interfaceBiologyGrid, &
               interfaceGrid, &
               verticalGrid, &
               nBioLayers, &
               nIceLayers, &
               nSnowLayers, &
               nAlgae, &
               nzAerosols, &
               nCategories, &
               nDOC, &
               nDIC, &
               nDON, &
               nDissolvedIron, &
               nParticulateIron, &
               basalIceMeltCategory(:,iCell), &
               surfaceIceMeltCategory(:,iCell), &
               congelationCategory(:,iCell), &
               snowiceFormationCategory(:,iCell), &
               seaSurfaceTemperature(iCell), &
               seaSurfaceSalinity(iCell), &
               snowfallRate(iCell), &
               snowMeltCategory(:,iCell), &
               oceanMixedLayerDepth(iCell), &
               initialSalinityProfile(:,iCell), &
               iceThicknessCategoryInitial(:,iCell), &
               oceanBioFluxes(:,iCell), &
               atmosBioFluxes(:,iCell), &
               iceAreaCategoryInitial(:,iCell), &
               iceVolumeCategoryInitial(:,iCell), &
               iceAreaCategory(1,:,iCell), &
               iceVolumeCategory(1,:,iCell), &
               snowVolumeCategory(1,:,iCell), &
               ciceTracerObject % tracerArrayCategory(:,:), &
               snowVolumeCategoryInitial(:,iCell), &
               config_use_skeletal_biochemistry, &
               maxAlgaeType, &
               nZBGCTracers, &
               abortFlag, &
               abortMessage)
          call colpkg_get_warnings(warnings)
          call column_write_warnings(warnings)

          ! code abort
          if (abortFlag) then
             COLUMN_ERROR_WRITE("column_biogeochemistry")
             COLUMN_ERROR_WRITE(trim(abortMessage))
             write(abortLocation,*) indexToCellID(iCell)
             call MPAS_dmpar_global_abort("MPAS-seaice: iCell "//trim(abortLocation)//": "//trim(abortMessage))
          endif

          call get_cice_tracer_array_category(block, ciceTracerObject, iCell, setGetPhysicsTracers, setGetBGCTracers)

          ! update
          totalSkeletalAlgae(iCell) = 0.0_RKIND

          do iCategory = 1, nCategories
             if (config_use_skeletal_biochemistry .and. iceAreaCategory(1,iCategory,iCell) > cicePuny) then
                do iAlgae = 1, nAlgae
                   totalSkeletalAlgae(iCell) = totalSkeletalAlgae(iCell) + &
                        skeletalAlgaeConc(iAlgae,iCategory,iCell) * &
                        iceAreaCategory(1,iCategory,iCell)
                enddo
             endif
          enddo ! iCategory

          do iCategory = 1, nCategories
             newlyFormedIce(iCategory,iCell) = 0
             if (newlyFormedIceLogical(iCategory)) newlyFormedIce(iCategory,iCell) = 1
          enddo ! iCategory

          if (.not. rayleighCriteria) rayleighCriteriaReal(iCell) = 0.0_RKIND

       enddo ! iCell

       deallocate(newlyFormedIceLogical)
       deallocate(oceanBioConcentrationsUsed)

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
    case ("gregorian_noleap")
       daysInYear = sum(daysInMonth)
    end select

  end subroutine get_days_in_year

!-----------------------------------------------------------------------
! Other routines
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_update_state
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 31st March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_column_update_state(domain, stateUpdateType, dt, iceAgeTimeOffset)

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
    call cice_column_aggregate(domain)
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

          COLUMN_ERROR_WRITE("cice_column_update_state: Unknown update type: "//trim(stateUpdateType))
          call MPAS_dmpar_global_abort("MPAS-seaice: cice_column_update_state: Unknown update type")

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

  end subroutine cice_column_update_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_aggregate
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 6th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_column_aggregate(domain)

    use ice_colpkg, only: colpkg_aggregate

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate, &
         icestate

    logical, pointer :: &
         config_use_column_biogeochemistry

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
         iCell

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    logical :: &
         setGetPhysicsTracers, &
         setGetBGCTracers

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)

       call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

       call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(tracers_aggregate, "iceVolumeCell", iceVolumeCell)
       call MPAS_pool_get_array(tracers_aggregate, "snowVolumeCell", snowVolumeCell)

       call MPAS_pool_get_array(icestate, "openWaterArea", openWaterArea)

       setGetPhysicsTracers = .true.
       setGetBGCTracers     = config_use_column_biogeochemistry

       do iCell = 1, nCellsSolve

          ! set the category tracer array
          call set_cice_tracer_array_category(block, ciceTracerObject, iCell, setGetPhysicsTracers, setGetBGCTracers)

          call colpkg_aggregate(&
               nCategories, &
               iceAreaCategory(1,:,iCell), &
               ciceTracerObject % tracerArrayCategory, & ! trcrn
               iceVolumeCategory(1,:,iCell), &
               snowVolumeCategory(1,:,iCell), &
               iceAreaCell(iCell), &
               ciceTracerObject % tracerArrayCell, & ! trcr
               iceVolumeCell(iCell), &
               snowVolumeCell(iCell), &
               openWaterArea(iCell), &
               ciceTracerObject % nTracers, &
               ciceTracerObject % parentIndex, & ! trcr_depend
               ciceTracerObject % firstAncestorMask, & ! trcr_base
               ciceTracerObject % ancestorNumber, & ! n_trcr_strata
               ciceTracerObject % ancestorIndices) ! nt_strata

          ! set the cell tracer array
          call get_cice_tracer_array_cell(block, ciceTracerObject, iCell, setGetPhysicsTracers, setGetBGCTracers)

       enddo ! iCell

       block => block % next
    end do

  end subroutine cice_column_aggregate

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_coupling_prep
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 6th April
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_column_coupling_prep(domain)

    use cice_constants, only: &
         cicePuny

    use ice_constants_colpkg, only: &
         rhofresh

    type(domain_type) :: domain

    type(block_type), pointer :: block

    logical, pointer :: &
         config_use_ocean_mixed_layer, &
         config_include_pond_freshwater_feedback, &
         config_use_column_biogeochemistry

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
         oceanDustIronFlux

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
         oceanDissolvedIronFlux

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory

    real(kind=RKIND), pointer :: &
         config_dt

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
         oceanBioFluxesAll

    call MPAS_pool_get_config(domain % configs, "config_use_ocean_mixed_layer", config_use_ocean_mixed_layer)
    call MPAS_pool_get_config(domain % configs, "config_dt", config_dt)
    call MPAS_pool_get_config(domain % configs, "config_include_pond_freshwater_feedback", config_include_pond_freshwater_feedback)
    call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

    if (config_use_ocean_mixed_layer) &
         call cice_column_ocean_mixed_layer(domain)

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
       call MPAS_pool_get_array(biogeochemistry, "oceanBioFluxes", oceanBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "oceanAlgaeFlux", oceanAlgaeFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDOCFlux", oceanDOCFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDICFlux", oceanDICFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDONFlux", oceanDONFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanParticulateIronFlux", oceanParticulateIronFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDissolvedIronFlux", oceanDissolvedIronFlux)

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
             if (solarZenithAngleCosine(iCell) > cicePuny) then

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
             pondFreshWaterFlux(iCell)  = pondFreshWaterFlux(iCell) * rhofresh / config_dt
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
           if (config_use_column_biogeochemistry) then

              oceanBioFluxesAll(:)              = 0.0_RKIND
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
              enddo

              ! Nitrate
              iBioData = iBioData+1
              oceanNitrateFlux(iCell) = oceanBioFluxesAll(iBioData)

              ! Polysaccharids and Lipids
              do iBioTracers = 1, maxDOCType
                  iBioData = iBioData+1
                  oceanDOCFlux(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
              enddo

              ! DIC
              do iBioTracers = 1, maxDICType
                  iBioData = iBioData+1
                  oceanDICFlux(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
              enddo

              ! Chlorophyll (not saved)
              iBioData = iBioData+maxAlgaeType

              ! Ammonium
              iBioData = iBioData+1
              oceanAmmoniumFlux(iCell) = oceanBioFluxesAll(iBioData)

              ! Silicate
              iBioData = iBioData+1
              oceanSilicateFlux(iCell) = oceanBioFluxesAll(iBioData)

              ! DMSPp
              iBioData = iBioData+1
              oceanDMSPpFlux(iCell) = oceanBioFluxesAll(iBioData)

              ! DMSPd
              iBioData = iBioData+1
              oceanDMSPdFlux(iCell) = oceanBioFluxesAll(iBioData)

              ! DMS
              iBioData = iBioData+1
              oceanDMSFlux(iCell) = oceanBioFluxesAll(iBioData)

              ! DON (Proteins)
              iBioData = iBioData+1
              do iBioTracers = 1, maxDONType
                  iBioData = iBioData+1
                  oceanDONFlux(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
              enddo

              ! Dissolved Iron
              do iBioTracers = 1, maxIronType
                  iBioData = iBioData+1
                  oceanDissolvedIronFlux(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
              enddo

              ! Particulate Iron
              do iBioTracers = 1, maxIronType
                  iBioData = iBioData+1
                  oceanParticulateIronFlux(iBioTracers,iCell) = oceanBioFluxesAll(iBioData)
              enddo

              ! Black Carbon (not saved)
              iBioData = iBioData + maxBCType

              ! Dust (combined)
              do iBioTracers = 1, maxDustType
                  iBioData = iBioData+1
                  oceanDustIronFlux(iCell) = oceanDustIronFlux(iCell) +  oceanBioFluxesAll(iBioData)
              enddo

              ! Humics
              iBioData = iBioData+1
              oceanHumicsFlux(iCell) = oceanBioFluxesAll(iBioData)

          endif ! config_use_column_biogeochemistry

       enddo ! iCell

       deallocate(oceanBioFluxesAll)

       block => block % next
    enddo

    call cice_column_scale_fluxes(domain)

  end subroutine cice_column_coupling_prep

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_scale_fluxes
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 6th April
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_column_scale_fluxes(domain)

    use cice_constants, only: &
         ciceStefanBoltzmann, &
         ciceFreshWaterFreezingPoint

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
         config_use_column_biogeochemistry

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

             if (config_use_column_biogeochemistry) then

                oceanNitrateFlux(iCell)        = oceanNitrateFlux(iCell)             * iceAreaInverse
                oceanSilicateFlux(iCell)       = oceanSilicateFlux(iCell)            * iceAreaInverse
                oceanAmmoniumFlux(iCell)       = oceanAmmoniumFlux(iCell)            * iceAreaInverse
                oceanDMSPpFlux(iCell)          = oceanDMSPpFlux(iCell)               * iceAreaInverse
                oceanDMSPdFlux(iCell)          = oceanDMSPdFlux(iCell)               * iceAreaInverse
                oceanDMSFlux(iCell)            = oceanDMSFlux(iCell)                 * iceAreaInverse
                oceanHumicsFlux(iCell)         = oceanHumicsFlux(iCell)              * iceAreaInverse
                oceanDustIronFlux(iCell)       = oceanDustIronFlux(iCell)            * iceAreaInverse

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
                  -ciceStefanBoltzmann * (seaFreezingTemperature(iCell) + ciceFreshWaterFreezingPoint)**4
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

             if (config_use_column_biogeochemistry) then

                oceanNitrateFlux(iCell)           = 0.0_RKIND
                oceanSilicateFlux(iCell)          = 0.0_RKIND
                oceanAmmoniumFlux(iCell)          = 0.0_RKIND
                oceanDMSPpFlux(iCell)             = 0.0_RKIND
                oceanDMSPdFlux(iCell)             = 0.0_RKIND
                oceanDMSFlux(iCell)               = 0.0_RKIND
                oceanHumicsFlux(iCell)            = 0.0_RKIND
                oceanDustIronFlux(iCell)          = 0.0_RKIND
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

  end subroutine cice_column_scale_fluxes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_ocean_mixed_layer
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 6th April
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_column_ocean_mixed_layer(domain)

    use ice_colpkg, only: &
         colpkg_atm_boundary, &
         colpkg_ocn_mixed_layer

    use ice_constants_colpkg, only: &
         albocn

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

          call colpkg_atm_boundary(&
               'ocn', &
               seaSurfaceTemperature(iCell), &
               airPotentialTemperature(iCell), &
               uAirVelocity(iCell), &
               vAirVelocity(iCell), &
               windSpeed(iCell), &
               airLevelHeight(iCell), &
               airSpecificHumidity(iCell), &
               airDensity(iCell), &
               airStressOceanU(iCell), &
               airStressOceanV(iCell), &
               atmosReferenceTemperature2mOcean(iCell), &
               atmosReferenceHumidity2mOcean(iCell), &
               potentialTemperatureDifference, &
               specificHumidityDifference, &
               latentTransferCoefficient, &
               sensibleTransferCoefficient, &
               airDragCoefficient(iCell), &
               airOceanDragCoefficientRatio(iCell))

          albedoVisibleDirectOcean(iCell)  = albocn
          albedoIRDirectOcean(iCell)       = albocn
          albedoVisibleDiffuseOcean(iCell) = albocn
          albedoIRDiffuseOcean(iCell)      = albocn

          call colpkg_ocn_mixed_layer(&
               albedoVisibleDirectOcean(iCell), &
               shortwaveVisibleDirectDown(iCell), &
               albedoIRDirectOcean(iCell), &
               shortwaveIRDirectDown(iCell), &
               albedoVisibleDiffuseOcean(iCell), &
               shortwaveVisibleDiffuseDown(iCell), &
               albedoIRDiffuseOcean(iCell), &
               shortwaveIRDiffuseDown(iCell), &
               seaSurfaceTemperature(iCell), &
               longwaveUpOcean(iCell), &
               sensibleHeatFluxOcean(iCell), &
               sensibleTransferCoefficient, &
               latentHeatFluxOcean(iCell), &
               latentTransferCoefficient, &
               evaporativeWaterFluxOcean(iCell), &
               longwaveDown(iCell), &
               potentialTemperatureDifference, &
               specificHumidityDifference, &
               iceAreaCell(iCell), &
               oceanHeatFlux(iCell), &
               oceanShortwaveFlux(iCell), &
               oceanMixedLayerDepth(iCell), &
               seaFreezingTemperature(iCell), &
               oceanHeatFluxConvergence(iCell), &
               freezingMeltingPotential(iCell), &
               config_dt)

       enddo ! iCell

       block => block % next
    enddo

    ! remove frazil from below ice shelves if were testing that
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

  end subroutine cice_column_ocean_mixed_layer

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
         config_use_column_biogeochemistry

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nCategories", nCategories)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nZBGCTracers", nZBGCTracers)

    call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

    ! get the number of CICE tracers in trcrn
    call init_column_tracer_object_tracer_number(domain, tracerObject)

    ! allocate the category tracer array
    allocate(tracerObject % tracerArrayCategory(tracerObject % nTracers,nCategories))

    ! allocate the cell tracer array
    allocate(tracerObject % tracerArrayCell(tracerObject % nTracers))

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
    if (config_use_column_biogeochemistry) then

       allocate(tracerObject % index_LayerIndexToDataArray(nZBGCTracers))
       allocate(tracerObject % index_LayerIndexToBioIndex(nZBGCTracers))

       ! set all indices for biogeochemistry including parent, ancestor and ancestor mask
       call init_column_tracer_object_for_biogeochemistry(domain, tracerObject)

    endif

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
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols, &
         config_use_brine, &
         config_use_column_biogeochemistry, &
         config_use_vertical_zsalinity, &
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
         config_use_zaerosols

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
    call MPAS_pool_get_config(domain % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_aerosols", config_use_aerosols)

    call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)
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
    if (config_use_cesm_meltponds .or. &
        config_use_level_meltponds .or. &
        config_use_topo_meltponds) &
         tracerObject % nTracers = tracerObject % nTracers + 2

    ! level or topo ponds
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) &
         tracerObject % nTracers = tracerObject % nTracers + 1

    ! aerosols
    if (config_use_aerosols) &
         tracerObject % nTracers = tracerObject % nTracers + nAerosols*4

    !-----------------------------------------------------------------------
    ! biogeochemistry
    !-----------------------------------------------------------------------

    if (config_use_column_biogeochemistry) then

       ! save tracer number without bio tracers counted
       tracerObject % nTracersNotBio = tracerObject % nTracers

       ! biogeochemical tracers
       tracerObject % nBioTracersLayer = 0

       ! brine height tracer
       if (config_use_brine) &
            tracerObject % nTracers = tracerObject % nTracers + 1

       ! vertical zSalinity
       if (config_use_vertical_zsalinity) then
          tracerObject % nTracers = tracerObject % nTracers + nBioLayers
          tracerObject % nBioTracersLayer = tracerObject % nBioTracersLayer + 1
       endif
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
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols
    integer :: &
         nTracers

    integer, pointer :: &
         nIceLayers, &
         nSnowLayers

    integer, parameter :: indexMissingValue = 0

    call MPAS_pool_get_config(domain % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(domain % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_aerosols", config_use_aerosols)

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

    if (config_use_cesm_meltponds .or. &
        config_use_level_meltponds .or. &
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
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols

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
    call MPAS_pool_get_config(domain % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_aerosols", config_use_aerosols)

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
       tracerObject % parentIndex(tracerObject % index_snowEnthalpy) = 2
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

    ! cesm melt ponds
    if (config_use_cesm_meltponds) then
       tracerObject % parentIndex(tracerObject % index_pondArea)  = 0
       tracerObject % parentIndex(tracerObject % index_pondDepth) = 2 + tracerObject % index_pondArea
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
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds

    call MPAS_pool_get_config(domain % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)

    ! initialize
    tracerObject % ancestorNumber = 0
    tracerObject % ancestorIndices = 0

    ! cesm meltponds
    if (config_use_cesm_meltponds) then

       ! melt pond depth
       tracerObject % ancestorNumber (tracerObject % index_pondDepth)   = 1
       tracerObject % ancestorIndices(tracerObject % index_pondDepth,1) = tracerObject % index_pondArea ! on melt pond area

    endif

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

  subroutine set_cice_tracer_array_category(block, tracerObject, iCell, setPhysicsTracers, setBGCTracers)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    integer, intent(in) :: &
         iCell

    logical, intent(in) :: &
         setPhysicsTracers, &
         setBGCTracers

    ! get physics tracers
    if (setPhysicsTracers) &
         call set_cice_physics_tracer_array_category(block, tracerObject, iCell)

    ! get BGC tracers
    if (setBGCTracers) &
         call set_cice_biogeochemistry_tracer_array_category(block, tracerObject, iCell)

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

  subroutine get_cice_tracer_array_category(block, tracerObject, iCell, getPhysicsTracers, getBGCTracers)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    integer, intent(in) :: &
         iCell

    logical, intent(in) :: &
         getPhysicsTracers, &
         getBGCTracers

    ! get physics tracers
    if (getPhysicsTracers) &
         call get_cice_physics_tracer_array_category(block, tracerObject, iCell)

    ! get BGC tracers
    if (getBGCTracers) &
         call get_cice_biogeochemistry_tracer_array_category(block, tracerObject, iCell)

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

  subroutine set_cice_tracer_array_cell(block, tracerObject, iCell, setPhysicsTracers, setBGCTracers)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    integer, intent(in) :: &
         iCell

    logical, intent(in) :: &
         setPhysicsTracers, &
         setBGCTracers

    ! get physics tracers
    if (setPhysicsTracers) &
         call set_cice_physics_tracer_array_cell(block, tracerObject, iCell)

    ! get BGC tracers
    if (setBGCTracers) &
         call set_cice_biogeochemistry_tracer_array_cell(block, tracerObject, iCell)

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

  subroutine get_cice_tracer_array_cell(block, tracerObject, iCell, getPhysicsTracers, getBGCTracers)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    integer, intent(in) :: &
         iCell

    logical, intent(in) :: &
         getPhysicsTracers, &
         getBGCTracers

    ! get physics tracers
    if (getPhysicsTracers) &
         call get_cice_physics_tracer_array_cell(block, tracerObject, iCell)

    ! get BGC tracers
    if (getBGCTracers) &
         call get_cice_biogeochemistry_tracer_array_cell(block, tracerObject, iCell)

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

  subroutine set_cice_physics_tracer_array_category(block, tracerObject, iCell)

    type(block_type), intent(in) :: &
         block

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols

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
         iceBodyAerosol

    integer :: &
         nTracers, &
         iAerosol

    call MPAS_pool_get_config(block % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(block % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(block % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(block % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)

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

    nTracers = 1

    ! surfaceTemperature
    tracerObject % tracerArrayCategory(nTracers,:) = surfaceTemperature(1,:,iCell)
    nTracers = nTracers + 1

    ! iceEnthalpy
    tracerObject % tracerArrayCategory(nTracers:nTracers+nIceLayers-1,:) = iceEnthalpy(:,:,iCell)
    nTracers = nTracers + nIceLayers

    ! snowEnthalpy
    tracerObject % tracerArrayCategory(nTracers:nTracers+nSnowLayers-1,:) = snowEnthalpy(:,:,iCell)
    nTracers = nTracers + nSnowLayers

    ! ice Salinity
    tracerObject % tracerArrayCategory(nTracers:nTracers+nIceLayers-1,:) = iceSalinity(:,:,iCell)
    nTracers = nTracers + nIceLayers

    ! iceAge
    if (config_use_ice_age) then
       tracerObject % tracerArrayCategory(nTracers,:) = iceAge(1,:,iCell)
       nTracers = nTracers + 1
    endif

    ! firstYearIceArea
    if (config_use_first_year_ice) then
       tracerObject % tracerArrayCategory(nTracers,:) = firstYearIceArea(1,:,iCell)
       nTracers = nTracers + 1
    endif

    ! level ice tracers
    if (config_use_level_ice) then
       tracerObject % tracerArrayCategory(nTracers,:) = levelIceArea(1,:,iCell)
       nTracers = nTracers + 1
       tracerObject % tracerArrayCategory(nTracers,:) = levelIceVolume(1,:,iCell)
       nTracers = nTracers + 1
    endif

    ! pond tracers
    if (config_use_cesm_meltponds .or. &
        config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       tracerObject % tracerArrayCategory(nTracers,:) = pondArea(1,:,iCell)
       nTracers = nTracers + 1
       tracerObject % tracerArrayCategory(nTracers,:) = pondDepth(1,:,iCell)
       nTracers = nTracers + 1
    endif

    ! level or topo ponds
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       tracerObject % tracerArrayCategory(nTracers,:) = pondLidThickness(1,:,iCell)
       nTracers = nTracers + 1
    end if

    ! aerosols
    if (config_use_aerosols) then
       do iAerosol = 1, nAerosols

          tracerObject % tracerArrayCategory(nTracers+4*(iAerosol-1)  ,:) = snowScatteringAerosol(iAerosol,:,iCell)
          tracerObject % tracerArrayCategory(nTracers+4*(iAerosol-1)+1,:) = snowBodyAerosol(iAerosol,:,iCell)
          tracerObject % tracerArrayCategory(nTracers+4*(iAerosol-1)+2,:) = iceScatteringAerosol(iAerosol,:,iCell)
          tracerObject % tracerArrayCategory(nTracers+4*(iAerosol-1)+3,:) = iceBodyAerosol(iAerosol,:,iCell)

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

  subroutine get_cice_physics_tracer_array_category(block, tracerObject, iCell)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols

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
         iceBodyAerosol

    integer :: &
         nTracers, &
         iAerosol

    call MPAS_pool_get_config(block % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(block % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(block % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(block % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)

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

    nTracers = 1

    ! surfaceTemperature
    surfaceTemperature(1,:,iCell) = tracerObject % tracerArrayCategory(nTracers,:)
    nTracers = nTracers + 1

    ! iceEnthalpy
    iceEnthalpy(:,:,iCell) = tracerObject % tracerArrayCategory(nTracers:nTracers+nIceLayers-1,:)
    nTracers = nTracers + nIceLayers

    ! snowEnthalpy
    snowEnthalpy(:,:,iCell) = tracerObject % tracerArrayCategory(nTracers:nTracers+nSnowLayers-1,:)
    nTracers = nTracers + nSnowLayers

    ! ice Salinity
    iceSalinity(:,:,iCell) = tracerObject % tracerArrayCategory(nTracers:nTracers+nIceLayers-1,:)
    nTracers = nTracers + nIceLayers

    ! iceAge
    if (config_use_ice_age) then
       iceAge(1,:,iCell) = tracerObject % tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
    endif

    ! firstYearIceArea
    if (config_use_first_year_ice) then
       firstYearIceArea(1,:,iCell) = tracerObject % tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
    endif

    ! level ice tracers
    if (config_use_level_ice) then
       levelIceArea(1,:,iCell) = tracerObject % tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
       levelIceVolume(1,:,iCell) = tracerObject % tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
    endif

    ! pond tracers
    if (config_use_cesm_meltponds .or. &
        config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       pondArea(1,:,iCell) = tracerObject % tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
       pondDepth(1,:,iCell) = tracerObject % tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
    endif

    ! level or topo ponds
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       pondLidThickness(1,:,iCell) = tracerObject % tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
    end if

    ! aerosols
    if (config_use_aerosols) then
       do iAerosol = 1, nAerosols

          snowScatteringAerosol(iAerosol,:,iCell) = tracerObject % tracerArrayCategory(nTracers+4*(iAerosol-1)  ,:)
          snowBodyAerosol(iAerosol,:,iCell)       = tracerObject % tracerArrayCategory(nTracers+4*(iAerosol-1)+1,:)
          iceScatteringAerosol(iAerosol,:,iCell)  = tracerObject % tracerArrayCategory(nTracers+4*(iAerosol-1)+2,:)
          iceBodyAerosol(iAerosol,:,iCell)        = tracerObject % tracerArrayCategory(nTracers+4*(iAerosol-1)+3,:)

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

  subroutine set_cice_physics_tracer_array_cell(block, tracerObject, iCell)

    type(block_type), intent(in) :: &
         block

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols

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
         iceBodyAerosolCell

    integer :: &
         nTracers, &
         iAerosol

    call MPAS_pool_get_config(block % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(block % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(block % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(block % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)

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

    nTracers = 1

    ! surfaceTemperature
    tracerObject % tracerArrayCell(nTracers) = surfaceTemperatureCell(iCell)
    nTracers = nTracers + 1

    ! iceEnthalpy
    tracerObject % tracerArrayCell(nTracers:nTracers+nIceLayers-1) = iceEnthalpyCell(:,iCell)
    nTracers = nTracers + nIceLayers

    ! snowEnthalpy
    tracerObject % tracerArrayCell(nTracers:nTracers+nSnowLayers-1) = snowEnthalpyCell(:,iCell)
    nTracers = nTracers + nSnowLayers

    ! ice Salinity
    tracerObject % tracerArrayCell(nTracers:nTracers+nIceLayers-1) = iceSalinityCell(:,iCell)
    nTracers = nTracers + nIceLayers

    ! iceAge
    if (config_use_ice_age) then
       tracerObject % tracerArrayCell(nTracers) = iceAgeCell(iCell)
       nTracers = nTracers + 1
    endif

    ! firstYearIceArea
    if (config_use_first_year_ice) then
       tracerObject % tracerArrayCell(nTracers) = firstYearIceAreaCell(iCell)
       nTracers = nTracers + 1
    endif

    ! level ice tracers
    if (config_use_level_ice) then
       tracerObject % tracerArrayCell(nTracers) = levelIceAreaCell(iCell)
       nTracers = nTracers + 1
       tracerObject % tracerArrayCell(nTracers) = levelIceVolumeCell(iCell)
       nTracers = nTracers + 1
    endif

    ! pond tracers
    if (config_use_cesm_meltponds .or. &
        config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       tracerObject % tracerArrayCell(nTracers) = pondAreaCell(iCell)
       nTracers = nTracers + 1
       tracerObject % tracerArrayCell(nTracers) = pondDepthCell(iCell)
       nTracers = nTracers + 1
    endif

    ! level or topo ponds
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       tracerObject % tracerArrayCell(nTracers) = pondLidThicknessCell(iCell)
       nTracers = nTracers + 1
    end if

    ! aerosols
    if (config_use_aerosols) then
       do iAerosol = 1, nAerosols

          tracerObject % tracerArrayCell(nTracers+4*(iAerosol-1)  ) = snowScatteringAerosolCell(iAerosol,iCell)
          tracerObject % tracerArrayCell(nTracers+4*(iAerosol-1)+1) = snowBodyAerosolCell(iAerosol,iCell)
          tracerObject % tracerArrayCell(nTracers+4*(iAerosol-1)+2) = iceScatteringAerosolCell(iAerosol,iCell)
          tracerObject % tracerArrayCell(nTracers+4*(iAerosol-1)+3) = iceBodyAerosolCell(iAerosol,iCell)

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

  subroutine get_cice_physics_tracer_array_cell(block, tracerObject, iCell)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols

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
         iceBodyAerosolCell

    integer :: &
         nTracers, &
         iAerosol

    call MPAS_pool_get_config(block % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(block % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(block % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(block % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)

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

    nTracers = 1

    ! surfaceTemperature
    surfaceTemperatureCell(iCell) = tracerObject % tracerArrayCell(nTracers)
    nTracers = nTracers + 1

    ! iceEnthalpy
    iceEnthalpyCell(:,iCell) = tracerObject % tracerArrayCell(nTracers:nTracers+nIceLayers-1)
    nTracers = nTracers + nIceLayers

    ! snowEnthalpy
    snowEnthalpyCell(:,iCell) = tracerObject % tracerArrayCell(nTracers:nTracers+nSnowLayers-1)
    nTracers = nTracers + nSnowLayers

    ! ice Salinity
    iceSalinityCell(:,iCell) = tracerObject % tracerArrayCell(nTracers:nTracers+nIceLayers-1)
    nTracers = nTracers + nIceLayers

    ! iceAge
    if (config_use_ice_age) then
       iceAgeCell(iCell) = tracerObject % tracerArrayCell(nTracers)
       nTracers = nTracers + 1
    endif

    ! firstYearIceArea
    if (config_use_first_year_ice) then
       firstYearIceAreaCell(iCell) = tracerObject % tracerArrayCell(nTracers)
       nTracers = nTracers + 1
    endif

    ! level ice tracers
    if (config_use_level_ice) then
       levelIceAreaCell(iCell) = tracerObject % tracerArrayCell(nTracers)
       nTracers = nTracers + 1
       levelIceVolumeCell(iCell) = tracerObject % tracerArrayCell(nTracers)
       nTracers = nTracers + 1
    endif

    ! pond tracers
    if (config_use_cesm_meltponds .or. &
        config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       pondAreaCell(iCell) = tracerObject % tracerArrayCell(nTracers)
       nTracers = nTracers + 1
       pondDepthCell(iCell) = tracerObject % tracerArrayCell(nTracers)
       nTracers = nTracers + 1
    endif

    ! level or topo ponds
    if (config_use_level_meltponds .or. &
        config_use_topo_meltponds) then
       pondLidThicknessCell(iCell) = tracerObject % tracerArrayCell(nTracers)
       nTracers = nTracers + 1
    end if

    ! aerosols
    if (config_use_aerosols) then
       do iAerosol = 1, nAerosols

          snowScatteringAerosolCell(iAerosol,iCell) = tracerObject % tracerArrayCell(nTracers+4*(iAerosol-1)  )
          snowBodyAerosolCell(iAerosol,iCell)       = tracerObject % tracerArrayCell(nTracers+4*(iAerosol-1)+1)
          iceScatteringAerosolCell(iAerosol,iCell)  = tracerObject % tracerArrayCell(nTracers+4*(iAerosol-1)+2)
          iceBodyAerosolCell(iAerosol,iCell)        = tracerObject % tracerArrayCell(nTracers+4*(iAerosol-1)+3)

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

  subroutine set_cice_biogeochemistry_tracer_array_category(block, tracerObject, iCell)

    type(block_type), intent(in) :: &
         block

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_skeletal_biochemistry, &
         config_use_vertical_biochemistry, &
         config_use_vertical_zsalinity, &
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
    call MPAS_pool_get_config(block % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)
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
         tracerObject % tracerArrayCategory(tracerObject % index_brineFraction,:) = brineFraction(1,:,iCell)

    if (config_use_skeletal_biochemistry) then

       ! algal nitrogen
       do iBioTracers = 1, nAlgae
          tracerObject % tracerArrayCategory(tracerObject % index_algaeConc(iBioTracers),:) = &
               skeletalAlgaeConc(iBioTracers,:,iCell)
       enddo

       ! nitrate
       if (config_use_nitrate) &
            tracerObject % tracerArrayCategory(tracerObject % index_nitrateConc,:) = skeletalNitrateConc(1,:,iCell)

       ! DOC
       if (config_use_carbon) then
          do iBioTracers = 1, nDOC
             tracerObject % tracerArrayCategory(tracerObject % index_DOCConc(iBioTracers),:) = skeletalDOCConc(iBioTracers,:,iCell)
          enddo

          ! DIC
          do iBioTracers = 1, nDIC
             tracerObject % tracerArrayCategory(tracerObject % index_DICConc(iBioTracers),:) = skeletalDICConc(iBioTracers,:,iCell)
          enddo
       endif

       ! DON
       if (config_use_DON) then
          do iBioTracers = 1, nDON
             tracerObject % tracerArrayCategory(tracerObject % index_DONConc(iBioTracers),:) = skeletalDONConc(iBioTracers,:,iCell)
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) &
            tracerObject % tracerArrayCategory(tracerObject % index_ammoniumConc,:) = skeletalAmmoniumConc(1,:,iCell)

       ! silicate
       if (config_use_silicate) &
            tracerObject % tracerArrayCategory(tracerObject % index_silicateConc,:) = skeletalSilicateConc(1,:,iCell)

       ! DMS, DMSPp, DMSPd
       if (config_use_DMS) then
          tracerObject % tracerArrayCategory(tracerObject % index_DMSConc,:) = skeletalDMSConc(1,:,iCell)
          tracerObject % tracerArrayCategory(tracerObject % index_DMSPpConc,:) = skeletalDMSPpConc(1,:,iCell)
          tracerObject % tracerArrayCategory(tracerObject % index_DMSPdConc,:) = skeletalDMSPdConc(1,:,iCell)
       endif

       ! nonreactive mobile tracer
       if (config_use_nonreactive) &
            tracerObject % tracerArrayCategory(tracerObject % index_nonreactiveConc,:) = skeletalNonreactiveConc(1,:,iCell)
       ! humic material
       if (config_use_humics) &
            tracerObject % tracerArrayCategory(tracerObject % index_humicsConc,:) = skeletalHumicsConc(1,:,iCell)

       ! Particulate and dissovled Iron
       if (config_use_iron) then
          do iBioTracers = 1, nParticulateIron
             tracerObject % tracerArrayCategory(tracerObject % index_particulateIronConc(iBioTracers),:) = &
                  skeletalParticulateIronConc(iBioTracers,:,iCell)
          enddo
          do iBioTracers = 1, nDissolvedIron
             tracerObject % tracerArrayCategory(tracerObject % index_dissolvedIronConc(iBioTracers),:) = &
                  skeletalDissolvedIronConc(iBioTracers,:,iCell)
          enddo
       endif

    elseif (config_use_vertical_tracers) then

       ! Fraction of biogeochemical tracer in the mobile phase
       do iLayers = 1, tracerObject % nBioTracers
          tracerObject % tracerArrayCategory(tracerObject % index_mobileFraction+iLayers-1,:) = mobileFraction(iLayers,:,iCell)
       enddo

       ! algal nitrogen
       if (config_use_vertical_biochemistry) then
          iBioCount = 0
          do iBioTracers = 1, nAlgae
             do iLayers = 1,nBioLayersP3
                iBiocount = iBiocount + 1
                tracerObject % tracerArrayCategory(tracerObject % index_algaeConc(iBioTracers)+iLayers-1,:) = &
                     verticalAlgaeConc(iBioCount,:,iCell)
             enddo
          enddo
       endif

       ! nitrate
       if (config_use_nitrate) then
          do iLayers = 1, nBioLayersP3
             tracerObject % tracerArrayCategory(tracerObject % index_nitrateConc + iLayers-1,:) = &
                  verticalNitrateConc(iLayers,:,iCell)
          enddo
       endif

       ! DOC
       if (config_use_carbon) then
          iBioCount = 0
          do iBioTracers = 1, nDOC
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerObject % tracerArrayCategory(tracerObject % index_DOCConc(iBioTracers) + iLayers-1,:) = &
                     verticalDOCConc(iBioCount,:,iCell)
             enddo
          enddo
          iBioCount = 0

          ! DIC
          do iBioTracers = 1, nDIC
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerObject % tracerArrayCategory(tracerObject % index_DICConc(iBioTracers) + iLayers-1,:) = &
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
                tracerObject % tracerArrayCategory(tracerObject % index_DONConc(iBioTracers) + iLayers-1,:) = &
                     verticalDONConc(iBioCount,:,iCell)
             enddo
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) then
          do iLayers = 1, nBioLayersP3
             tracerObject % tracerArrayCategory(tracerObject % index_ammoniumConc + iLayers-1,:) = &
                  verticalAmmoniumConc(iLayers,:,iCell)
          enddo
       endif

       ! silicate
       if (config_use_silicate) then
          do iLayers = 1, nBioLayersP3
             tracerObject % tracerArrayCategory(tracerObject % index_silicateConc+iLayers-1,:) = &
                  verticalSilicateConc(iLayers,:,iCell)
          enddo
       endif

       ! DMS, DMSPp, DMSPd
       if (config_use_DMS) then
          do iLayers = 1, nBioLayersP3
             tracerObject % tracerArrayCategory(tracerObject % index_DMSConc+iLayers-1,:) = verticalDMSConc(iLayers,:,iCell)
             tracerObject % tracerArrayCategory(tracerObject % index_DMSPpConc+iLayers-1,:) = verticalDMSPpConc(iLayers,:,iCell)
             tracerObject % tracerArrayCategory(tracerObject % index_DMSPdConc+iLayers-1,:) = verticalDMSPdConc(iLayers,:,iCell)
          enddo
       endif

       ! nonreactive purely mobile tracers
       if (config_use_nonreactive) then
          do iLayers = 1, nBioLayersP3
             tracerObject % tracerArrayCategory(tracerObject % index_nonreactiveConc+iLayers-1,:) = &
                  verticalNonreactiveConc(iLayers,:,iCell)
          enddo
       endif

       ! humic material
       if (config_use_humics) then
          do iLayers = 1, nBioLayersP3
             tracerObject % tracerArrayCategory(tracerObject % index_humicsConc+iLayers-1,:) = verticalHumicsConc(iLayers,:,iCell)
          enddo
       endif

       ! particulate and dissolved Iron
       if (config_use_iron) then
          iBioCount = 0
          do iBioTracers = 1, nParticulateIron
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerObject % tracerArrayCategory(tracerObject % index_particulateIronConc(iBioTracers)+iLayers-1,:) = &
                     verticalParticulateIronConc(iBioCount,:,iCell)
             enddo
          enddo
          iBioCount = 0
          do iBioTracers = 1, nDissolvedIron
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerObject % tracerArrayCategory(tracerObject % index_dissolvedIronConc(iBioTracers)+iLayers-1,:) = &
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
                tracerObject % tracerArrayCategory(tracerObject % index_verticalAerosolsConc(iBioTracers)+iLayers-1,:) = &
                     verticalAerosolsConc(iBioCount,:,iCell)
             enddo
          enddo
       endif

       ! salinity used with BL99 thermodynamics
       if (config_use_vertical_zsalinity) then
          do iLayers = 1, nBioLayers
             tracerObject % tracerArrayCategory(tracerObject % index_verticalSalinity+iLayers-1,:) = &
                  verticalSalinity(iLayers,:,iCell)
          enddo
       endif
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

  subroutine get_cice_biogeochemistry_tracer_array_category(block, tracerObject, iCell)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_skeletal_biochemistry, &
         config_use_vertical_biochemistry, &
         config_use_vertical_zsalinity, &
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
    call MPAS_pool_get_config(block % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)
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
         brineFraction(1,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_brineFraction,:)

    if (config_use_skeletal_biochemistry) then

       ! algal nitrogen
       do iBioTracers = 1, nAlgae
          skeletalAlgaeConc(iBioTracers,:,iCell) = &
               tracerObject % tracerArrayCategory(tracerObject % index_algaeConc(iBioTracers),:)
       enddo

       ! nitrate
       if (config_use_nitrate) &
            skeletalNitrateConc(1,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_nitrateConc,:)

       if (config_use_carbon) then

          ! DOC
          do iBioTracers = 1, nDOC
             skeletalDOCConc(iBioTracers,:,iCell) = &
                  tracerObject % tracerArrayCategory(tracerObject % index_DOCConc(iBioTracers),:)
          enddo

          ! DIC
          do iBioTracers = 1, nDIC
             skeletalDICConc(iBioTracers,:,iCell) = &
                  tracerObject % tracerArrayCategory(tracerObject % index_DICConc(iBioTracers),:)
          enddo
       endif

       ! DON
       if (config_use_DON) then
          do iBioTracers = 1, nDON
             skeletalDONConc(iBioTracers,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_DONConc(iBioTracers),:)
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) &
            skeletalAmmoniumConc(1,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_ammoniumConc,:)
       ! silicate
       if (config_use_silicate) &
            skeletalSilicateConc(1,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_silicateConc,:)
       ! DNS, DMSPp, DMSPd
       if (config_use_DMS) then
          skeletalDMSConc(1,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_DMSConc,:)
          skeletalDMSPpConc(1,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_DMSPpConc,:)
          skeletalDMSPdConc(1,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_DMSPdConc,:)
       endif

       ! nonreactive tracer
       if (config_use_nonreactive) &
            skeletalNonreactiveConc(1,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_nonreactiveConc,:)
       ! humic material
       if (config_use_humics) &
            skeletalHumicsConc(1,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_humicsConc,:)

       if (config_use_iron) then

          ! Particulate Iron
          do iBioTracers = 1, nParticulateIron
             skeletalParticulateIronConc(iBioTracers,:,iCell) = &
                  tracerObject % tracerArrayCategory(tracerObject % index_particulateIronConc(iBioTracers),:)
          enddo

          ! Dissolved Iron
          do iBioTracers = 1, nDissolvedIron
             skeletalDissolvedIronConc(iBioTracers,:,iCell) = &
                  tracerObject % tracerArrayCategory(tracerObject % index_dissolvedIronConc(iBioTracers),:)
          enddo
       endif

    elseif (config_use_vertical_tracers) then

       ! fraction of biogeochemical tracer in the mobile phase
       do iLayers = 1, tracerObject % nBioTracers
          mobileFraction(iLayers,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_mobileFraction+iLayers-1,:)
       enddo

       if (config_use_vertical_biochemistry) then
          iBioCount = 0

          ! algal nitrogen
          do iBioTracers = 1, nAlgae
             do iLayers = 1,nBioLayersP3
                iBiocount = iBiocount + 1
                verticalAlgaeConc(iBioCount,:,iCell) = &
                     tracerObject % tracerArrayCategory(tracerObject % index_algaeConc(iBioTracers)+iLayers-1,:)
             enddo
          enddo
       endif

       ! nitrate
       if (config_use_nitrate) then
          do iLayers = 1, nBioLayersP3
             verticalNitrateConc(iLayers,:,iCell) = &
                  tracerObject % tracerArrayCategory(tracerObject % index_nitrateConc + iLayers-1,:)
          enddo
       endif

       if (config_use_carbon) then
          iBioCount = 0

          ! DOC
          do iBioTracers = 1, nDOC
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDOCConc(iBioCount,:,iCell) = &
                     tracerObject % tracerArrayCategory(tracerObject % index_DOCConc(iBioTracers) + iLayers-1,:)
             enddo
          enddo
          iBioCount = 0

          ! DIC
          do iBioTracers = 1, nDIC
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDICConc(iBioCount,:,iCell) = &
                     tracerObject % tracerArrayCategory(tracerObject % index_DICConc(iBioTracers) + iLayers-1,:)
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
                     tracerObject % tracerArrayCategory(tracerObject % index_DONConc(iBioTracers) + iLayers-1,:)
             enddo
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) then
          do iLayers = 1, nBioLayersP3
             verticalAmmoniumConc(iLayers,:,iCell) = &
                  tracerObject % tracerArrayCategory(tracerObject % index_ammoniumConc + iLayers-1,:)
          enddo
       endif

       ! silicate
       if (config_use_silicate) then
          do iLayers = 1, nBioLayersP3
             verticalSilicateConc(iLayers,:,iCell) = &
                  tracerObject % tracerArrayCategory(tracerObject % index_silicateConc+iLayers-1,:)
          enddo
       endif

       ! DMS, DMSPp, DMSPd
       if (config_use_DMS) then
          do iLayers = 1, nBioLayersP3
             verticalDMSConc(iLayers,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_DMSConc+iLayers-1,:)
             verticalDMSPpConc(iLayers,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_DMSPpConc+iLayers-1,:)
             verticalDMSPdConc(iLayers,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_DMSPdConc+iLayers-1,:)
          enddo
       endif

       ! nonreactive tracer
       if (config_use_nonreactive) then
          do iLayers = 1, nBioLayersP3
             verticalNonreactiveConc(iLayers,:,iCell) = &
                  tracerObject % tracerArrayCategory(tracerObject % index_nonreactiveConc+iLayers-1,:)
          enddo
       endif

       ! humic material
       if (config_use_humics) then
          do iLayers = 1, nBioLayersP3
             verticalHumicsConc(iLayers,:,iCell) = tracerObject % tracerArrayCategory(tracerObject % index_humicsConc+iLayers-1,:)
          enddo
       endif

       if (config_use_iron) then
          iBioCount = 0

          ! particulate iron
          do iBioTracers = 1, nParticulateIron
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalParticulateIronConc(iBioCount,:,iCell) = &
                     tracerObject % tracerArrayCategory(tracerObject % index_particulateIronConc(iBioTracers)+iLayers-1,:)
             enddo
          enddo
          iBioCount = 0

          ! dissolved iron
          do iBioTracers = 1, nDissolvedIron
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDissolvedIronConc(iBioCount,:,iCell) = &
                     tracerObject % tracerArrayCategory(tracerObject % index_dissolvedIronConc(iBioTracers)+iLayers-1,:)
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
                     tracerObject % tracerArrayCategory(tracerObject % index_verticalAerosolsConc(iBioTracers)+iLayers-1,:)
             enddo
          enddo
       endif

       ! salinity used with BL99 thermodynamics
       if (config_use_vertical_zsalinity) then
          do iLayers = 1, nBioLayers
             verticalSalinity(iLayers,:,iCell) = &
                  tracerObject % tracerArrayCategory(tracerObject % index_verticalSalinity+iLayers-1,:)
          enddo
       endif
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

  subroutine set_cice_biogeochemistry_tracer_array_cell(block, tracerObject, iCell)

    type(block_type), intent(in) :: &
         block

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_skeletal_biochemistry, &
         config_use_vertical_biochemistry, &
         config_use_vertical_zsalinity, &
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
         verticalSalinityCell

    integer :: &
         iBioTracers, &
         iBioCount, &
         iLayers

    call MPAS_pool_get_config(block % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)
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
    call MPAS_pool_get_array(tracers_aggregate, "verticalSalinityCell", verticalSalinityCell)
    call MPAS_pool_get_array(tracers_aggregate, "brineFractionCell", brineFractionCell)

    ! biogeochemistry
    ! brine height fraction
    if (config_use_brine) &
         tracerObject % tracerArrayCell(tracerObject % index_brineFraction)  = brineFractionCell(iCell)

    if (config_use_skeletal_biochemistry) then

       ! algal nitrogen
       do iBioTracers = 1, nAlgae
          tracerObject % tracerArrayCell(tracerObject % index_algaeConc(iBioTracers)) = skeletalAlgaeConcCell(iBioTracers,iCell)
       enddo

       ! nitrate
       if (config_use_nitrate) &
            tracerObject % tracerArrayCell(tracerObject % index_nitrateConc)  = skeletalNitrateConcCell(1,iCell)

       if (config_use_carbon) then

          ! DOC
          do iBioTracers = 1, nDOC
             tracerObject % tracerArrayCell(tracerObject % index_DOCConc(iBioTracers)) = skeletalDOCConcCell(iBioTracers,iCell)
          enddo

          ! DIC
          do iBioTracers = 1, nDIC
             tracerObject % tracerArrayCell(tracerObject % index_DICConc(iBioTracers))  = skeletalDICConcCell(iBioTracers,iCell)
          enddo
       endif

       ! DON
       if (config_use_DON) then
          do iBioTracers = 1, nDON
             tracerObject % tracerArrayCell(tracerObject % index_DONConc(iBioTracers))  = skeletalDONConcCell(iBioTracers,iCell)
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) &
            tracerObject % tracerArrayCell(tracerObject % index_ammoniumConc)  = skeletalAmmoniumConcCell(1,iCell)

       ! silicate
       if (config_use_silicate) &
            tracerObject % tracerArrayCell(tracerObject % index_silicateConc)  = skeletalSilicateConcCell(1,iCell)

       ! DMS, DMSPp, DMSPd
       if (config_use_DMS) then
          tracerObject % tracerArrayCell(tracerObject % index_DMSConc)  = skeletalDMSConcCell(1,iCell)
          tracerObject % tracerArrayCell(tracerObject % index_DMSPpConc) = skeletalDMSPpConcCell(1,iCell)
          tracerObject % tracerArrayCell(tracerObject % index_DMSPdConc) = skeletalDMSPdConcCell(1,iCell)
       endif

       ! nonreactive tracer
       if (config_use_nonreactive) &
            tracerObject % tracerArrayCell(tracerObject % index_nonreactiveConc) = skeletalNonreactiveConcCell(1,iCell)
       ! humic material
       if (config_use_humics) &
            tracerObject % tracerArrayCell(tracerObject % index_humicsConc)  = skeletalHumicsConcCell(1,iCell)

       if (config_use_iron) then

          ! particulate iron
          do iBioTracers = 1, nParticulateIron
             tracerObject % tracerArrayCell(tracerObject % index_particulateIronConc(iBioTracers)) = &
                  skeletalParticulateIronConcCell(iBioTracers,iCell)
          enddo

          ! dissolved iron
          do iBioTracers = 1, nDissolvedIron
             tracerObject % tracerArrayCell(tracerObject % index_dissolvedIronConc(iBioTracers)) = &
                  skeletalDissolvedIronConcCell(iBioTracers,iCell)
          enddo
       endif

    elseif (config_use_vertical_tracers) then

       if (config_use_vertical_biochemistry) then
          iBioCount = 0

          ! algal nitrogen
          do iBioTracers = 1, nAlgae
             do iLayers = 1,nBioLayersP3
                iBiocount = iBiocount + 1
                tracerObject % tracerArrayCell(tracerObject % index_algaeConc(iBioTracers)+iLayers-1) = &
                     verticalAlgaeConcCell(iBioCount,iCell)
             enddo
          enddo
       endif

       ! nitrate
       if (config_use_nitrate) then
          do iLayers = 1, nBioLayersP3
             tracerObject % tracerArrayCell(tracerObject % index_nitrateConc + iLayers-1) = verticalNitrateConcCell(iLayers,iCell)
          enddo
       endif

       if (config_use_carbon) then
          iBioCount = 0

          ! DOC
          do iBioTracers = 1, nDOC
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerObject % tracerArrayCell(tracerObject % index_DOCConc(iBioTracers) + iLayers-1) = &
                     verticalDOCConcCell(iBioCount,iCell)
             enddo
          enddo
          iBioCount = 0

          ! DIC
          do iBioTracers = 1, nDIC
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerObject % tracerArrayCell(tracerObject % index_DICConc(iBioTracers) + iLayers-1) = &
                     verticalDICConcCell(iBioCount,iCell)
             enddo
          enddo
       endif

       ! DON
       if (config_use_DON) then
          iBioCount = 0
          do iBioTracers = 1, nDON
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerObject % tracerArrayCell(tracerObject % index_DONConc(iBioTracers) + iLayers-1) = &
                     verticalDONConcCell(iBioCount,iCell)
             enddo
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) then
          do iLayers = 1, nBioLayersP3
             tracerObject % tracerArrayCell(tracerObject % index_ammoniumConc + iLayers-1) = &
                  verticalAmmoniumConcCell(iLayers,iCell)
          enddo
       endif

       ! silicate
       if (config_use_silicate) then
          do iLayers = 1, nBioLayersP3
             tracerObject % tracerArrayCell(tracerObject % index_silicateConc+iLayers-1) = verticalSilicateConcCell(iLayers,iCell)
          enddo
       endif

       ! DMS, DMSPp, DMSPd
       if (config_use_DMS) then
          do iLayers = 1, nBioLayersP3
             tracerObject % tracerArrayCell(tracerObject % index_DMSConc+iLayers-1)  = verticalDMSConcCell(iLayers,iCell)
             tracerObject % tracerArrayCell(tracerObject % index_DMSPpConc+iLayers-1)  = verticalDMSPpConcCell(iLayers,iCell)
             tracerObject % tracerArrayCell(tracerObject % index_DMSPdConc+iLayers-1)  = verticalDMSPdConcCell(iLayers,iCell)
          enddo
       endif

       ! nonreactive
       if (config_use_nonreactive) then
          do iLayers = 1, nBioLayersP3
             tracerObject % tracerArrayCell(tracerObject % index_nonreactiveConc+iLayers-1) = &
                  verticalNonreactiveConcCell(iLayers,iCell)
          enddo
       endif

       ! humic material
       if (config_use_humics) then
          do iLayers = 1, nBioLayersP3
             tracerObject % tracerArrayCell(tracerObject % index_humicsConc+iLayers-1)  = verticalHumicsConcCell(iLayers,iCell)
          enddo
       endif

       if (config_use_iron) then
          iBioCount = 0

          ! particulate iron
          do iBioTracers = 1, nParticulateIron
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerObject % tracerArrayCell(tracerObject % index_particulateIronConc(iBioTracers)+iLayers-1) = &
                     verticalParticulateIronConcCell(iBioCount,iCell)
             enddo
          enddo
          iBioCount = 0

          ! dissolved iron
          do iBioTracers = 1, nDissolvedIron
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerObject % tracerArrayCell(tracerObject % index_dissolvedIronConc(iBioTracers)+iLayers-1) = &
                     verticalDissolvedIronConcCell(iBioCount,iCell)
             enddo
          enddo
       endif

       ! black carbon and dust aerosols
       if (config_use_zaerosols) then
          iBioCount = 0
          do iBioTracers = 1, nzAerosols
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                tracerObject % tracerArrayCell(tracerObject % index_verticalAerosolsConc(iBioTracers)+iLayers-1) = &
                     verticalAerosolsConcCell(iBioCount,iCell)
             enddo
          enddo
       endif

       ! salinity for use with BL99 thermodynamics
       if (config_use_vertical_zsalinity) then
          do iLayers = 1, nBioLayers
             tracerObject % tracerArrayCell(tracerObject % index_verticalSalinity+iLayers-1)  = verticalSalinityCell(iLayers,iCell)
          enddo
       endif
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

  subroutine get_cice_biogeochemistry_tracer_array_cell(block, tracerObject, iCell)

    type(block_type), intent(inout) :: &
         block

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    integer, intent(in) :: &
         iCell

    logical, pointer :: &
         config_use_skeletal_biochemistry, &
         config_use_vertical_biochemistry, &
         config_use_vertical_zsalinity, &
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
         verticalSalinityCell

    integer :: &
         iBioTracers, &
         iBioCount, &
         iLayers

    call MPAS_pool_get_config(block % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(block % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)
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
    call MPAS_pool_get_array(tracers_aggregate, "verticalSalinityCell", verticalSalinityCell)
    call MPAS_pool_get_array(tracers_aggregate, "brineFractionCell", brineFractionCell)

    ! biogeochemistry
    ! brine height fraction
    if (config_use_brine) &
         brineFractionCell(iCell) = tracerObject % tracerArrayCell(tracerObject % index_brineFraction)

    if (config_use_skeletal_biochemistry) then

       ! algal nitrogen
       do iBioTracers = 1, nAlgae
          skeletalAlgaeConcCell(iBioTracers,iCell) = tracerObject % tracerArrayCell(tracerObject % index_algaeConc(iBioTracers))
       enddo

       ! nitrate
       if (config_use_nitrate) &
            skeletalNitrateConcCell(1,iCell) = tracerObject % tracerArrayCell(tracerObject % index_nitrateConc)

       if (config_use_carbon) then

          ! DOC
          do iBioTracers = 1, nDOC
             skeletalDOCConcCell(iBioTracers,iCell) =  tracerObject % tracerArrayCell(tracerObject % index_DOCConc(iBioTracers))
          enddo

          ! DIC
          do iBioTracers = 1, nDIC
             skeletalDICConcCell(iBioTracers,iCell) = tracerObject % tracerArrayCell(tracerObject % index_DICConc(iBioTracers))
          enddo
       endif

       ! DON
       if (config_use_DON) then
          do iBioTracers = 1, nDON
             skeletalDONConcCell(iBioTracers,iCell) = tracerObject % tracerArrayCell(tracerObject % index_DONConc(iBioTracers))
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) &
            skeletalAmmoniumConcCell(1,iCell) = tracerObject % tracerArrayCell(tracerObject % index_ammoniumConc)

       ! silicate
       if (config_use_silicate) &
            skeletalSilicateConcCell(1,iCell) = tracerObject % tracerArrayCell(tracerObject % index_silicateConc)

       ! DMS
       if (config_use_DMS) then
          skeletalDMSConcCell(1,iCell) = tracerObject % tracerArrayCell(tracerObject % index_DMSConc)
          skeletalDMSPpConcCell(1,iCell) = tracerObject % tracerArrayCell(tracerObject % index_DMSPpConc)
          skeletalDMSPdConcCell(1,iCell) = tracerObject % tracerArrayCell(tracerObject % index_DMSPdConc)
       endif

       ! nonreactive tracer
       if (config_use_nonreactive) &
            skeletalNonreactiveConcCell(1,iCell) = tracerObject % tracerArrayCell(tracerObject % index_nonreactiveConc)
       ! humic material
       if (config_use_humics) &
            skeletalHumicsConcCell(1,iCell) = tracerObject % tracerArrayCell(tracerObject % index_humicsConc)

       if (config_use_iron) then

          ! particulate iron
          do iBioTracers = 1, nParticulateIron
             skeletalParticulateIronConcCell(iBioTracers,iCell) = &
                  tracerObject % tracerArrayCell(tracerObject % index_particulateIronConc(iBioTracers))
          enddo

          ! dissolved iron
          do iBioTracers = 1, nDissolvedIron
             skeletalDissolvedIronConcCell(iBioTracers,iCell) = &
                  tracerObject % tracerArrayCell(tracerObject % index_dissolvedIronConc(iBioTracers))
          enddo
       endif

    elseif (config_use_vertical_tracers) then

       if (config_use_vertical_biochemistry) then
          iBioCount = 0

          ! algal nitrogen
          do iBioTracers = 1, nAlgae
             do iLayers = 1,nBioLayersP3
                iBiocount = iBiocount + 1
                verticalAlgaeConcCell(iBioCount,iCell) = &
                     tracerObject % tracerArrayCell(tracerObject % index_algaeConc(iBioTracers)+iLayers-1)
             enddo
          enddo
       endif

       ! nitrate
       if (config_use_nitrate) then
          do iLayers = 1, nBioLayersP3
             verticalNitrateConcCell(iLayers,iCell) = tracerObject % tracerArrayCell(tracerObject % index_nitrateConc + iLayers-1)
          enddo
       endif

       if (config_use_carbon) then
          iBioCount = 0

          ! DOC
          do iBioTracers = 1, nDOC
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDOCConcCell(iBioCount,iCell) = &
                     tracerObject % tracerArrayCell(tracerObject % index_DOCConc(iBioTracers) + iLayers-1)
             enddo
          enddo
          iBioCount = 0

          ! DIC
          do iBioTracers = 1, nDIC
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDICConcCell(iBioCount,iCell) = &
                     tracerObject % tracerArrayCell(tracerObject % index_DICConc(iBioTracers) + iLayers-1)
             enddo
          enddo
       endif

       ! DON
       if (config_use_DON) then
          iBioCount = 0
          do iBioTracers = 1, nDON
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDONConcCell(iBioCount,iCell) = &
                     tracerObject % tracerArrayCell(tracerObject % index_DONConc(iBioTracers) + iLayers-1)
             enddo
          enddo
       endif

       ! ammonium
       if (config_use_ammonium) then
          do iLayers = 1, nBioLayersP3
             verticalAmmoniumConcCell(iLayers,iCell) = &
                  tracerObject % tracerArrayCell(tracerObject % index_ammoniumConc + iLayers-1)
          enddo
       endif

       ! silicate
       if (config_use_silicate) then
          do iLayers = 1, nBioLayersP3
             verticalSilicateConcCell(iLayers,iCell) = tracerObject % tracerArrayCell(tracerObject % index_silicateConc+iLayers-1)
          enddo
       endif

       ! DMS, DMSPp, DMSPd
       if (config_use_DMS) then
          do iLayers = 1, nBioLayersP3
             verticalDMSConcCell(iLayers,iCell) = tracerObject % tracerArrayCell(tracerObject % index_DMSConc+iLayers-1)
             verticalDMSPpConcCell(iLayers,iCell) = tracerObject % tracerArrayCell(tracerObject % index_DMSPpConc+iLayers-1)
             verticalDMSPdConcCell(iLayers,iCell) = tracerObject % tracerArrayCell(tracerObject % index_DMSPdConc+iLayers-1)
          enddo
       endif

       ! nonreactive tracer
       if (config_use_nonreactive) then
          do iLayers = 1, nBioLayersP3
             verticalNonreactiveConcCell(iLayers,iCell) = &
                  tracerObject % tracerArrayCell(tracerObject % index_nonreactiveConc+iLayers-1)
          enddo
       endif

       ! humic material
       if (config_use_humics) then
          do iLayers = 1, nBioLayersP3
             verticalHumicsConcCell(iLayers,iCell) = tracerObject % tracerArrayCell(tracerObject % index_humicsConc+iLayers-1)
          enddo
       endif

       if (config_use_iron) then
          iBioCount = 0

          ! particulate iron
          do iBioTracers = 1, nParticulateIron
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalParticulateIronConcCell(iBioCount,iCell) = &
                     tracerObject % tracerArrayCell(tracerObject % index_particulateIronConc(iBioTracers)+iLayers-1)
             enddo
          enddo
          iBioCount = 0

          ! dissolved iron
          do iBioTracers = 1, nDissolvedIron
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalDissolvedIronConcCell(iBioCount,iCell) = &
                     tracerObject % tracerArrayCell(tracerObject % index_dissolvedIronConc(iBioTracers)+iLayers-1)
             enddo
          enddo
       endif

       ! black carbon and dust aerosols
       if (config_use_zaerosols) then
          iBioCount = 0
          do iBioTracers = 1, nzAerosols
             do iLayers = 1,nBioLayersP3
                iBioCount = iBioCount + 1
                verticalAerosolsConcCell(iBioCount,iCell) = &
                     tracerObject % tracerArrayCell(tracerObject % index_verticalAerosolsConc(iBioTracers)+iLayers-1)
             enddo
          enddo
       endif

       ! salinity for use with BL99 thermodynamics
       if (config_use_vertical_zsalinity) then
          do iLayers = 1, nBioLayers
             verticalSalinityCell(iLayers,iCell) = tracerObject % tracerArrayCell(tracerObject % index_verticalSalinity+iLayers-1)
          enddo
       endif
    endif

  end subroutine get_cice_biogeochemistry_tracer_array_cell

!-----------------------------------------------------------------------
! Init CICE parameters
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_package_parameters
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2nd Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_package_parameters(domain, tracerObject)

    type(domain_type), intent(inout) :: domain

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    ! check column configs
    call check_column_package_configs(domain)

    ! set the tracer flags
    call init_column_package_tracer_flags(domain)

    ! set the tracer numbers
    call init_column_package_tracer_numbers(tracerObject)

    ! set the tracers indices
    call init_column_package_tracer_indices(tracerObject)

    ! set the column parameters
    call init_column_package_configs(domain)

  end subroutine init_column_package_parameters

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  check_column_package_configs
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 5th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine check_column_package_configs(domain)

    use cice_constants, only: &
         cicePuny

    type(domain_type), intent(inout) :: &
         domain

    integer, pointer :: &
         nCategories

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
         config_ocean_heat_transfer_type, &
         config_sea_freezing_temperature_type

    logical, pointer :: &
         config_calc_surface_stresses, &
         config_calc_surface_temperature, &
         config_use_form_drag, &
         config_use_level_ice, &
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_vertical_zsalinity, &
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
         config_use_column_biogeochemistry

    logical :: &
         use_meltponds

    integer :: &
         nPondSchemesActive

    real(kind=RKIND), pointer :: &
         config_max_meltwater_retained_fraction, &
         config_min_meltwater_retained_fraction, &
         config_snow_to_ice_transition_depth

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nCategories", nCategories)

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
    call MPAS_pool_get_config(domain % configs, "config_calc_surface_stresses", config_calc_surface_stresses)
    call MPAS_pool_get_config(domain % configs, "config_calc_surface_temperature", config_calc_surface_temperature)
    call MPAS_pool_get_config(domain % configs, "config_max_meltwater_retained_fraction", config_max_meltwater_retained_fraction)
    call MPAS_pool_get_config(domain % configs, "config_min_meltwater_retained_fraction", config_min_meltwater_retained_fraction)
    call MPAS_pool_get_config(domain % configs, "config_snow_to_ice_transition_depth", config_snow_to_ice_transition_depth)
    call MPAS_pool_get_config(domain % configs, "config_use_form_drag", config_use_form_drag)
    call MPAS_pool_get_config(domain % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_ocean_heat_transfer_type", config_ocean_heat_transfer_type)
    call MPAS_pool_get_config(domain % configs, "config_sea_freezing_temperature_type", config_sea_freezing_temperature_type)
    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)
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

    !-----------------------------------------------------------------------
    ! Check values
    !-----------------------------------------------------------------------

    ! check config_thermodynamics_type value
    if (.not. (trim(config_thermodynamics_type) == "zero layer" .or. &
               trim(config_thermodynamics_type) == "BL99" .or. &
               trim(config_thermodynamics_type) == "mushy")) then
       call config_error("config_thermodynamics_type", config_thermodynamics_type, "'zero layer', 'BL99' or 'mushy'")
    endif

    ! check config_heat_conductivity_type value
    if (.not. (trim(config_heat_conductivity_type) == "MU71" .or. &
               trim(config_heat_conductivity_type) == "bubbly")) then
       call config_error("config_heat_conductivity_type", config_heat_conductivity_type, "'MU71' or 'bubbly'")
    endif

    ! check config_shortwave_type value
    if (.not. (trim(config_shortwave_type) == "ccsm3" .or. &
               trim(config_shortwave_type) == "dEdd")) then
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
    if (.not. (trim(config_atmos_boundary_method) == "ccsm3" .or. &
               trim(config_atmos_boundary_method) == "constant")) then
       call config_error("config_atmos_boundary_method", config_atmos_boundary_method, "'ccsm3' or 'constant'")
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

    !-----------------------------------------------------------------------
    ! Check combinations
    !-----------------------------------------------------------------------

    ! check only single meltpond option on
    nPondSchemesActive = 0
    if (config_use_cesm_meltponds)  nPondSchemesActive = nPondSchemesActive + 1
    if (config_use_level_meltponds) nPondSchemesActive = nPondSchemesActive + 1
    if (config_use_topo_meltponds)  nPondSchemesActive = nPondSchemesActive + 1
    if (nPondSchemesActive > 1) then
       COLUMN_ERROR_WRITE('More than one melt pond scheme active')
       call MPAS_dmpar_global_abort('MPAS-seaice: More than one melt pond scheme active')
    endif

    ! check for itd remapping with only one category
    if (nCategories == 1 .and. trim(config_itd_conversion_type) == "linear remap") then
       COLUMN_ERROR_WRITE('Remapping the ITD is not allowed for nCategories=1')
       COLUMN_ERROR_WRITE("Use config_itd_conversion_type = 'delta function' with config_category_bounds_type = 'original'")
       COLUMN_ERROR_WRITE("or for column configurations use config_category_bounds_type = 'single category'")
       call MPAS_dmpar_global_abort('MPAS-seaice: Incompatible nCategories and config_itd_conversion_type')
    endif

    ! check itd and category bounds discrepancy
    if (nCategories /= 1 .and. trim(config_category_bounds_type) == 'single category') then
       COLUMN_ERROR_WRITE("nCategories /= 1 .and. config_category_bounds_type = 'single category'")
       call MPAS_dmpar_global_abort('MPAS-seaice: Incompatible nCategories and config_category_bounds_type')
    endif

    ! check config_snow_to_ice_transition_depth and level ponds
    if (config_use_level_meltponds .and. abs(config_snow_to_ice_transition_depth) > cicePuny) then
       COLUMN_ERROR_WRITE("config_use_level_meltponds = .true. and config_snow_to_ice_transition_depth /= 0")
       call MPAS_dmpar_global_abort('MPAS-seaice: Incompatible config_use_level_meltponds and config_snow_to_ice_transition_depth')
    endif

    ! check cesm ponds and freezing lids inconsistency
    if (config_use_cesm_meltponds .and. trim(config_pond_refreezing_type) /= "cesm") then
       COLUMN_ERROR_WRITE("config_use_cesm_meltponds = .true. and config_pond_refreezing_type /= 'cesm'")
       call MPAS_dmpar_global_abort('MPAS-seaice: Incompatible config_use_cesm_meltponds and config_pond_refreezing_type')
    endif

    ! check dEdd shortwave if using ponds
    use_meltponds = (config_use_cesm_meltponds .or. config_use_level_meltponds .or. config_use_topo_meltponds)
    if (trim(config_shortwave_type) /= 'dEdd' .and. use_meltponds .and. config_calc_surface_temperature) then
       COLUMN_ERROR_WRITE("config_shortwave_type) /= 'dEdd' .and. use_meltponds = .true. .")
       COLUMN_ERROR_WRITE("and. config_calc_surface_temperature ==.true.")
       call MPAS_dmpar_global_abort(&
            'MPAS-seaice: Incompatible config_shortwave_type, use_meltponds and config_calc_surface_temperature')
    endif

    ! check range of config_min_meltwater_retained_fraction and config_max_meltwater_retained_fraction
    if (config_min_meltwater_retained_fraction < 0.0_RKIND .or. &
        config_min_meltwater_retained_fraction > 1.0_RKIND) then
       COLUMN_ERROR_WRITE('config_min_meltwater_retained_fraction out of bounds')
       call MPAS_dmpar_global_abort('MPAS-seaice: config_min_meltwater_retained_fraction out of bounds')
    endif
    if (config_max_meltwater_retained_fraction < 0.0_RKIND .or. &
        config_max_meltwater_retained_fraction > 1.0_RKIND) then
       COLUMN_ERROR_WRITE('config_max_meltwater_retained_fraction out of bounds')
       call MPAS_dmpar_global_abort('MPAS-seaice: config_max_meltwater_retained_fraction out of bounds')
    endif

    ! check not mushy physics and dont calculate surface temperature
    if (trim(config_thermodynamics_type) == "mushy" .and. .not. config_calc_surface_temperature) then
       COLUMN_ERROR_WRITE("config_thermodynamics_type = 'mushy' and config_calc_surface_temperature = .false.")
       call MPAS_dmpar_global_abort('MPAS-seaice: Incompatible config_thermodynamics_type and config_calc_surface_temperature')
    endif

    ! check not form drag with constant atmosphere boundary method
    if (config_use_form_drag .and. trim(config_atmos_boundary_method) == "constant") then
       COLUMN_ERROR_WRITE("config_use_form_drag = .true. and config_atmos_boundary_method = 'constant'")
       call MPAS_dmpar_global_abort('MPAS-seaice: Incompatible config_use_form_drag and config_atmos_boundary_method')
    endif

    ! check not form drag with not calculating surface stresses
    if (config_use_form_drag .and. .not. config_calc_surface_stresses) then
       COLUMN_ERROR_WRITE("config_use_form_drag = .true. and config_calc_surface_stresses = .false.")
       call MPAS_dmpar_global_abort('MPAS-seaice: Incompatible config_use_form_drag and config_calc_surface_stresses')
    endif

    ! check am not using form drag with cesm ponds
    if (config_use_form_drag .and. config_use_cesm_meltponds) then
       COLUMN_ERROR_WRITE("config_use_form_drag = .true. and config_use_cesm_meltponds = .true.")
       call MPAS_dmpar_global_abort('MPAS-seaice: Incompatible config_use_form_drag and config_use_cesm_meltponds')
    endif

    ! check using form drag but not level ice
    if (config_use_form_drag .and. .not. config_use_level_ice) then
       COLUMN_ERROR_WRITE("config_use_form_drag = .true. and config_use_level_ice = .false.")
       call MPAS_dmpar_global_abort('MPAS-seaice: Incompatible config_use_form_drag and config_use_level_ice')
    endif

    ! check form drag and ocean heat flux type
    if (.not. config_use_form_drag .and. trim(config_ocean_heat_transfer_type) == "Cdn_ocn") then
       COLUMN_ERROR_WRITE("config_use_form_drag = .false. and config_ocean_heat_transfer_type == 'Cdn_ocn'")
       call MPAS_dmpar_global_abort('MPAS-seaice: Incompatible config_use_form_drag and config_ocean_heat_transfer_type')
    endif

    ! check thermodynamic type and sea freezing temperature type
    if (trim(config_thermodynamics_type) == "BL99" .and. trim(config_sea_freezing_temperature_type) /= "linear_salt") then
       COLUMN_ERROR_WRITE("config_thermodynamics_type == 'BL99' and config_sea_freezing_temperature_type /= 'linear_salt'")
       call MPAS_dmpar_global_abort(&
            'MPAS-seaice: Incompatible config_thermodynamics_type and config_sea_freezing_temperature_type')
    endif
    if (trim(config_thermodynamics_type) == "mushy" .and. trim(config_sea_freezing_temperature_type) /= "mushy") then
       COLUMN_ERROR_WRITE("config_thermodynamics_type == 'mushy' and config_sea_freezing_temperature_type /= 'mushy'")
       call MPAS_dmpar_global_abort(&
            'MPAS-seaice: Incompatible config_thermodynamics_type and config_sea_freezing_temperature_type')
    endif

    ! check biogeochemistry flags:
    if (.not. config_use_column_biogeochemistry .and. (config_use_brine .or. config_use_vertical_zsalinity .or. &
        config_use_vertical_biochemistry .or. config_use_shortwave_bioabsorption .or. config_use_vertical_tracers .or. &
        config_use_skeletal_biochemistry .or. config_use_nitrate .or. config_use_carbon .or. config_use_chlorophyll .or. &
        config_use_ammonium .or. config_use_silicate .or. config_use_DMS .or. config_use_nonreactive .or. config_use_humics .or. &
        config_use_DON .or. config_use_iron .or. config_use_modal_aerosols .or. config_use_zaerosols)) then
       COLUMN_ERROR_WRITE("config_use_column_biogeochemistry = false. All biogeochemistry namelist flags must also be false")
       call MPAS_dmpar_global_abort &
            ('MPAS-seaice: Incompatible config_use_column_biogeochemistry and true biogeochemistry flag')
    endif

    ! check vertical zSalinity requirements
    if (config_use_vertical_zsalinity .and. ((.not. config_use_brine) .or. &
         (.not. (trim(config_thermodynamics_type) == "BL99")))) then
       COLUMN_ERROR_WRITE("vertical zSalinity requires config_use_brine = true and 'BL99' ")
       call MPAS_dmpar_global_abort &
            ('MPAS-seaice: Incompatible config_use_vertical_zsalinity, config_use_brine, and config_thermodynamics_type')
    endif

    ! check that vertical bio tracers use brine height
    if ((config_use_vertical_biochemistry .or. config_use_zaerosols) .and. &
         (.not. config_use_brine .or. .not. config_use_vertical_tracers )) then
       COLUMN_ERROR_WRITE("vertical biochemistry and zaerosols require config_use_brine and config_use_vertical_tracer = true")
       call MPAS_dmpar_global_abort &
            ('MPAS-seaice: Incompatible config_use_vertical_biochemistry, '// &
             'config_use_zaerosols, config_use_brine, and config_use_vertical_tracer')
    endif

    ! check that vertical bio tracers and not used with skeletal bio tracers
    if (config_use_vertical_tracers  .and.  config_use_skeletal_biochemistry) then
       COLUMN_ERROR_WRITE("vertical bio tracers and skeletal bio tracers cannot both be true")
       call MPAS_dmpar_global_abort('MPAS-seaice: Incompatible config_use_vertical_tracers, config_use_skeletal_biochemistry')
    endif

    ! check that the shortwave scheme and bioabsorption is consistent
    if (config_use_shortwave_bioabsorption .and. .not. (trim(config_shortwave_type) == "dEdd")) then
       COLUMN_ERROR_WRITE("shortwave bioabsorption requires config_shortwave_type ==  'dEdd'")
       call MPAS_dmpar_global_abort('MPAS-seaice: Incompatible config_use_shortwave, config_shortwave_type')
    endif

    ! check that nitrate is true for biogeochemistry
    if ((config_use_vertical_biochemistry .or. config_use_skeletal_biochemistry) .and. .not. config_use_nitrate) then
       COLUMN_ERROR_WRITE("biochemistry needs at the very least config_use_nitrate = true")
       call MPAS_dmpar_global_abort&
            ('MPAS-seaice: Incompatible config_use_vertical_biochemistry, '// &
             'config_use_skeletal_biochemistry, and config_use_nitrate')
    endif

  end subroutine check_column_package_configs

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_package_tracer_flags
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2nd Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_package_tracer_flags(domain)

    !use ice_colpkg_tracers, only: &
    !     tr_iage      , & ! if .true., use age tracer
    !     tr_FY        , & ! if .true., use first-year area tracer
    !     tr_lvl       , & ! if .true., use level ice tracer
    !     tr_pond      , & ! if .true., use melt pond tracer
    !     tr_pond_cesm , & ! if .true., use cesm pond tracer
    !     tr_pond_lvl  , & ! if .true., use level-ice pond tracer
    !     tr_pond_topo , & ! if .true., use explicit topography-based ponds
    !     tr_aero      , & ! if .true., use aerosol tracers
    !     tr_brine         ! if .true., brine height differs from ice thickness

    use ice_colpkg, only: &
         colpkg_init_tracer_flags

    type(domain_type), intent(inout) :: domain

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols, &
         config_use_brine, &
         config_use_vertical_zsalinity, &
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
         config_use_skeletal_biochemistry

    logical :: &
         use_meltponds, &
         use_nitrogen

    call MPAS_pool_get_config(domain % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(domain % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_aerosols", config_use_aerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)
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

    use_nitrogen = .false.
    if (config_use_skeletal_biochemistry .or. config_use_vertical_biochemistry) &
         use_nitrogen = .true.

    use_meltponds = (config_use_cesm_meltponds .or. config_use_level_meltponds .or. config_use_topo_meltponds)

    call colpkg_init_tracer_flags(&
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         use_meltponds, &
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols, &
         config_use_brine, &
         config_use_vertical_zsalinity, &
         config_use_zaerosols, &
         config_use_nitrate, &
         use_nitrogen, &
         config_use_DON, &
         config_use_carbon, &
         config_use_chlorophyll, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_iron, &
         config_use_humics, &
         config_use_nonreactive)

    !tr_iage      = config_use_ice_age
    !tr_FY        = config_use_first_year_ice
    !tr_lvl       = config_use_level_ice
    !tr_pond      = use_meltponds
    !tr_pond_cesm = config_use_cesm_meltponds
    !tr_pond_lvl  = config_use_level_meltponds
    !tr_pond_topo = config_use_topo_meltponds
    !tr_aero      = config_use_aerosols
    !tr_brine     = config_use_brine
    !tr_bgc_S     = config_use_vertical_zsalinity
    !tr_zaero     = config_use_zaerosols
    !tr_bgc_Nit   = config_use_nitrate
    !tr_bgc_N     = use_nitrogen
    !tr_bgc_DON   = config_use_DON
    !tr_bgc_C     = config_use_carbon
    !tr_bgc_chl   = config_use_chlorophyll
    !tr_bgc_Am    = config_use_ammonium
    !tr_bgc_Sil   = config_use_silicate
    !tr_bgc_DMS   = config_use_DMS
    !tr_bgc_Fe    = config_use_iron
    !tr_bgc_hum   = config_use_humics
    !tr_bgc_PON   = config_use_nonreactive

  end subroutine init_column_package_tracer_flags

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_package_tracer_numbers
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 9th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_package_tracer_numbers(tracerObject)

    !use ice_colpkg_tracers, only: &
    !     ntrcr, &
    !     nbtrcr, &
    !     nbtrcr_sw

    use ice_colpkg, only: &
         colpkg_init_tracer_numbers

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    call colpkg_init_tracer_numbers(&
         tracerObject % nTracers, &
         tracerObject % nBioTracers, &
         tracerObject % nBioTracersShortwave)

    !ntrcr     = tracerObject % nTracers
    !nbtrcr    = tracerObject % nBioTracers
    !nbtrcr_sw = tracerObject % nBioTracersShortwave

  end subroutine init_column_package_tracer_numbers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_package_tracer_indices
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 5th Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_package_tracer_indices(tracerObject)

    !use ice_colpkg_tracers, only: &
    !     nt_Tsfc, &       ! ice/snow temperature
    !     nt_qice, &       ! volume-weighted ice enthalpy (in layers)
    !     nt_qsno, &       ! volume-weighted snow enthalpy (in layers)
    !     nt_sice, &       ! volume-weighted ice bulk salinity (CICE grid layers)
    !     nt_fbri, &       ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
    !     nt_iage, &       ! volume-weighted ice age
    !     nt_FY, &         ! area-weighted first-year ice area
    !     nt_alvl, &       ! level ice area fraction
    !     nt_vlvl, &       ! level ice volume fraction
    !     nt_apnd, &       ! melt pond area fraction
    !     nt_hpnd, &       ! melt pond depth
    !     nt_ipnd, &       ! melt pond refrozen lid thickness
    !     nt_aero, &       ! starting index for aerosols in ice
    !     nt_fbri, &       ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
    !     nt_bgc_Nit, &    ! nutrients
    !     nt_bgc_Am, &     !
    !     nt_bgc_Sil, &    !
    !     nt_bgc_DMSPp, &  ! trace gases (skeletal layer)
    !     nt_bgc_DMSPd, &  !
    !     nt_bgc_DMS, &    !
    !     nt_bgc_PON, &    ! zooplankton and detritus
    !     nt_bgc_hum, &    ! humic material
    !                      ! bio layer indicess
    !     nlt_bgc_Nit, &   ! nutrients
    !     nlt_bgc_Am, &    !
    !     nlt_bgc_Sil, &   !
    !     nlt_bgc_DMSPp, & ! trace gases (skeletal layer)
    !     nlt_bgc_DMSPd, & !
    !     nlt_bgc_DMS, &   !
    !     nlt_bgc_PON, &   ! zooplankton and detritus
    !     nlt_bgc_hum, &   ! humic material
    !     nlt_chl_sw, &    ! points to total chla in trcrn_sw
    !     nt_zbgc_frac, &  ! fraction of tracer in the mobile phase
    !     nt_bgc_S, &      ! Bulk salinity in fraction ice with dynamic salinity (Bio grid)
    !     nt_bgc_N, &      ! diatoms, phaeocystis, pico/small
    !     nt_bgc_C, &      ! diatoms, phaeocystis, pico/small
    !     nt_bgc_chl, &    ! diatoms, phaeocystis, pico/small
    !     nlt_bgc_N, &     ! diatoms, phaeocystis, pico/small
    !     nlt_bgc_C, &     ! diatoms, phaeocystis, pico/small
    !     nlt_bgc_chl, &   ! diatoms, phaeocystis, pico/small
    !     nt_bgc_DOC, &    ! dissolved organic carbon
    !     nlt_bgc_DOC, &   ! dissolved organic carbon
    !     nt_bgc_DON, &    ! dissolved organic nitrogen
    !     nlt_bgc_DON, &   ! dissolved organic nitrogen
    !     nt_bgc_DIC, &    ! dissolved inorganic carbon
    !     nlt_bgc_DIC, &   ! dissolved inorganic carbon
    !     nt_bgc_Fed, &    ! dissolved iron
    !     nt_bgc_Fep, &    ! particulate iron
    !     nlt_bgc_Fed, &   ! dissolved iron
    !     nlt_bgc_Fep, &   ! particulate iron
    !     nt_zaero, &      ! black carbon and other aerosols
    !     nlt_zaero, &     ! black carbon and other aerosols
    !     nlt_zaero_sw     ! black carbon and other aerosols

    use ice_colpkg, only: &
         colpkg_init_tracer_indices

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    call colpkg_init_tracer_indices(&
         tracerObject % index_surfaceTemperature, &
         tracerObject % index_iceEnthalpy, &
         tracerObject % index_snowEnthalpy, &
         tracerObject % index_iceSalinity, &
         tracerObject % index_brineFraction, &
         tracerObject % index_iceAge, &
         tracerObject % index_firstYearIceArea, &
         tracerObject % index_levelIceArea, &
         tracerObject % index_levelIceVolume, &
         tracerObject % index_pondArea, &
         tracerObject % index_pondDepth, &
         tracerObject % index_pondLidThickness, &
         tracerObject % index_aerosols, &
         tracerObject % index_verticalAerosolsConc, &
         tracerObject % index_algaeConc, &
         tracerObject % index_algalCarbon, &
         tracerObject % index_algalChlorophyll, &
         tracerObject % index_DOCConc, &
         tracerObject % index_DONConc, &
         tracerObject % index_DICConc, &
         tracerObject % index_dissolvedIronConc, &
         tracerObject % index_particulateIronConc, &
         tracerObject % index_nitrateConc, &
         tracerObject % index_ammoniumConc, &
         tracerObject % index_silicateConc, &
         tracerObject % index_DMSPpConc, &
         tracerObject % index_DMSPdConc, &
         tracerObject % index_DMSConc, &
         tracerObject % index_humicsConc, &
         tracerObject % index_nonreactiveConc, &
         tracerObject % index_verticalAerosolsConcLayer, &
         tracerObject % index_algaeConcLayer, &
         tracerObject % index_algalCarbonLayer, &
         tracerObject % index_algalChlorophyllLayer, &
         tracerObject % index_DOCConcLayer, &
         tracerObject % index_DONConcLayer, &
         tracerObject % index_DICConcLayer, &
         tracerObject % index_dissolvedIronConcLayer, &
         tracerObject % index_particulateIronConcLayer, &
         tracerObject % index_nitrateConcLayer, &
         tracerObject % index_ammoniumConcLayer, &
         tracerObject % index_silicateConcLayer, &
         tracerObject % index_DMSPpConcLayer, &
         tracerObject % index_DMSPdConcLayer, &
         tracerObject % index_DMSConcLayer, &
         tracerObject % index_humicsConcLayer, &
         tracerObject % index_nonreactiveConcLayer, &
         tracerObject % index_mobileFraction, &
         tracerObject % index_verticalSalinity, &
         tracerObject % index_chlorophyllShortwave, &
         tracerObject % index_verticalAerosolsConcShortwave, &
         tracerObject % nAlgaeIndex, &
         tracerObject % nAlgalCarbonIndex, &
         tracerObject % nAlgalChlorophyllIndex, &
         tracerObject % nDOCIndex, &
         tracerObject % nDONIndex, &
         tracerObject % nDICIndex, &
         tracerObject % nDissolvedIronIndex, &
         tracerObject % nParticulateIronIndex, &
         tracerObject % nzAerosolsIndex, &
         tracerObject % index_LayerIndexToDataArray, &
         tracerObject % index_LayerIndexToBioIndex, &
         tracerObject % nBioTracers)

    !nt_Tsfc       = tracerObject % index_surfaceTemperature
    !nt_qice       = tracerObject % index_iceEnthalpy
    !nt_qsno       = tracerObject % index_snowEnthalpy
    !nt_sice       = tracerObject % index_iceSalinity
    !nt_iage       = tracerObject % index_iceAge
    !nt_FY         = tracerObject % index_firstYearIceArea
    !nt_alvl       = tracerObject % index_levelIceArea
    !nt_vlvl       = tracerObject % index_levelIceVolume
    !nt_apnd       = tracerObject % index_pondArea
    !nt_hpnd       = tracerObject % index_pondDepth
    !nt_ipnd       = tracerObject % index_pondLidThickness
    !nt_aero       = tracerObject % index_aerosols
    !nt_fbri       = tracerObject % index_brineFraction
    !nt_zaeros     = tracerObject % index_verticalAerosolsConc
    !nt_bgc_N      = tracerObject % index_algaeConc
    !nt_bgc_C      = tracerObject % index_algalCarbon
    !nt_bgc_chl    = tracerObject % index_algalChlorophyll
    !nt_bgc_DOC    = tracerObject % index_DOCConc
    !nt_bgc_DON    = tracerObject % index_DONConc
    !nt_bgc_DIC    = tracerObject % index_DICConc
    !nt_bgc_Fed    = tracerObject % index_dissolvedIronConc
    !nt_bgc_Fep    = tracerObject % index_particulateIronConc
    !nt_bgc_Nit    = tracerObject % index_nitrateConc
    !nt_bgc_Am     = tracerObject % index_ammoniumConc
    !nt_bgc_Sil    = tracerObject % index_silicateConc
    !nt_bgc_DMSPp  = tracerObject % index_DMSPpConc
    !nt_bgc_DMSPd  = tracerObject % index_DMSPdConc
    !nt_bgc_DMS    = tracerObject % index_DMSConc
    !nt_bgc_hum    = tracerObject % index_humicsConc
    !nt_bgc_PON    = tracerObject % index_nonreactiveConc
    !nlt_zaero     = tracerObject % index_verticalAerosolsConcLayer
    !nlt_bgc_N     = tracerObject % index_algaeConcLayer
    !nlt_bgc_C     = tracerObject % index_algalCarbonLayer
    !nlt_bgc_chl   = tracerObject % index_algalChlorophyllLayer
    !nlt_bgc_DOC   = tracerObject % index_DOCConcLayer
    !nlt_bgc_DON   = tracerObject % index_DONConcLayer
    !nlt_bgc_DIC   = tracerObject % index_DICConcLayer
    !nlt_bgc_Fed   = tracerObject % index_dissolvedIronConcLayer
    !nlt_bgc_Fep   = tracerObject % index_particulateIronConcLayer
    !nlt_bgc_Nit   = tracerObject % index_nitrateConcLayer
    !nlt_bgc_Am    = tracerObject % index_ammoniumConcLayer
    !nlt_bgc_Sil   = tracerObject % index_silicateConcLayer
    !nlt_bgc_DMSPp = tracerObject % index_DMSPpConcLayer
    !nlt_bgc_DMSPd = tracerObject % index_DMSPdConcLayer
    !nlt_bgc_DMS   = tracerObject % index_DMSConcLayer
    !nlt_bgc_hum   = tracerObject % index_humicsConcLayer
    !nlt_bgc_PON   = tracerObject % index_nonreactiveConcLayer
    !nt_zbgc_frac  = tracerObject % index_mobileFraction
    !nt_zbgc_S     = tracerObject % index_verticalSalinity
    !nlt_chl_sw    = tracerObject % index_chlorophyllShortwave
    !nlt_zaero_sw  = tracerObject % index_verticalAerosolsConcShortwave
    !max_algae     = tracerObject % nAlgaeIndex
    !max_algae     = tracerObject % nAlgalCarbonIndex
    !max_algae     = tracerObject % nAlgalChlorophyllIndex
    !max_doc       = tracerObject % nDOCIndex
    !max_don       = tracerObject % nDONIndex
    !max_dic       = tracerObject % nDICIndex
    !max_fe        = tracerObject % nDissolvedIronIndex
    !max_fe        = tracerObject % nParticulateIronIndex
    !max_aero      = tracerObject % nzAerosolsIndex
    !bio_index_o   = tracerObject % index_LayerIndexToDataArray
    !bio_index     = tracerObject % index_LayerIndexToBioIndex
    !nbtrcr        = tracerObject % nBioTracers

  end subroutine init_column_package_tracer_indices

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_package_configs
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2nd Feburary 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_package_configs(domain)

    !use ice_colpkg_shared, only: &
    !     ktherm, &
    !     conduct, &
    !     fbot_xfer_type, &
    !     heat_capacity, &
    !     calc_Tsfc, &
    !     ustar_min, &
    !     a_rapid_mode, &
    !     Rac_rapid_mode, &
    !     aspect_rapid_mode, &
    !     dSdt_slow_mode, &
    !     phi_c_slow_mode, &
    !     phi_i_mushy, &
    !     shortwave, &
    !     albedo_type, &
    !     albicev, &
    !     albicei, &
    !     albsnowv, &
    !     albsnowi, &
    !     ahmax, &
    !     R_ice, &
    !     R_pnd, &
    !     R_snw, &
    !     dT_mlt, &
    !     rsnw_mlt, &
    !     kalg, &
    !     kstrength, &
    !     krdg_partic, &
    !     krdg_redist, &
    !     mu_rdg, &
    !     Cf, &
    !     atmbndy, &
    !     calc_strair, &
    !     formdrag, &
    !     highfreq, &
    !     natmiter, &
    !     oceanmixed_ice, &
    !     tfrz_option, &
    !     kitd, &
    !     kcatbound, &
    !     hs0, &
    !     frzpnd, &
    !     dpscale, &
    !     rfracmin, &
    !     rfracmax, &
    !     pndaspect, &
    !     hs1, &
    !     hp1
    !     bgc_flux_type, &
    !     z_tracers, &
    !     scale_bgc, &
    !     solve_zbgc, &
    !     dEdd_algae, &
    !     modal_aero, &
    !     skl_bgc, &
    !     solve_zsal, &
    !     grid_o, &
    !     l_sk, &
    !     grid_o_t, &
    !     initbio_frac, &
    !     frazil_scav, &
    !     grid_oS, &
    !     l_skS, &
    !     phi_snow, &
    !     ratio_Si2N_diatoms, &
    !     ratio_Si2N_sp     , &
    !     ratio_Si2N_phaeo  , &
    !     ratio_S2N_diatoms , &
    !     ratio_S2N_sp      , &
    !     ratio_S2N_phaeo   , &
    !     ratio_Fe2C_diatoms, &
    !     ratio_Fe2C_sp     , &
    !     ratio_Fe2C_phaeo  , &
    !     ratio_Fe2N_diatoms, &
    !     ratio_Fe2N_sp     , &
    !     ratio_Fe2N_phaeo  , &
    !     ratio_Fe2DON      , &
    !     ratio_Fe2DOC_s    , &
    !     ratio_Fe2DOC_l    , &
    !     fr_resp           , &
    !     tau_min           , &
    !     tau_max           , &
    !     algal_vel         , &
    !     R_dFe2dust        , &
    !     dustFe_sol        , &
    !     chlabs_diatoms    , &
    !     chlabs_sp         , &
    !     chlabs_phaeo      , &
    !     alpha2max_low_diatoms , &
    !     alpha2max_low_sp      , &
    !     alpha2max_low_phaeo   , &
    !     beta2max_diatoms , &
    !     beta2max_sp      , &
    !     beta2max_phaeo   , &
    !     mu_max_diatoms   , &
    !     mu_max_sp        , &
    !     mu_max_phaeo     , &
    !     grow_Tdep_diatoms, &
    !     grow_Tdep_sp     , &
    !     grow_Tdep_phaeo  , &
    !     fr_graze_diatoms , &
    !     fr_graze_sp      , &
    !     fr_graze_phaeo   , &
    !     mort_pre_diatoms , &
    !     mort_pre_sp      , &
    !     mort_pre_phaeo   , &
    !     mort_Tdep_diatoms, &
    !     mort_Tdep_sp     , &
    !     mort_Tdep_phaeo  , &
    !     k_exude_diatoms  , &
    !     k_exude_sp       , &
    !     k_exude_phaeo    , &
    !     K_Nit_diatoms    , &
    !     K_Nit_sp         , &
    !     K_Nit_phaeo      , &
    !     K_Am_diatoms     , &
    !     K_Am_sp          , &
    !     K_Am_phaeo       , &
    !     K_Sil_diatoms    , &
    !     K_Sil_sp         , &
    !     K_Sil_phaeo      , &
    !     K_Fe_diatoms     , &
    !     K_Fe_sp          , &
    !     K_Fe_phaeo       , &
    !     f_don_protein    , &
    !     kn_bac_protein   , &
    !     f_don_Am_protein , &
    !     f_doc_s            , &
    !     f_doc_l            , &
    !     f_exude_s          , &
    !     f_exude_l          , &
    !     k_bac_s            , &
    !     k_bac_l            , &
    !     T_max              , &
    !     fsal               , &
    !     op_dep_min         , &
    !     fr_graze_s         , &
    !     fr_graze_e         , &
    !     fr_mort2min        , &
    !     fr_dFe             , &
    !     k_nitrif           , &
    !     t_iron_conv        , &
    !     max_loss           , &
    !     max_dfe_doc1       , &
    !     fr_resp_s          , &
    !     y_sk_DMS           , &
    !     t_sk_conv          , &
    !     t_sk_ox            , &
    !     algaltype_diatoms  , &
    !     algaltype_sp       , &
    !     algaltype_phaeo    , &
    !     nitratetype        , &
    !     ammoniumtype       , &
    !     silicatetype       , &
    !     dmspptype          , &
    !     dmspdtype          , &
    !     humtype            , &
    !     doctype_s          , &
    !     doctype_l          , &
    !     dontype_protein    , &
    !     fedtype_1          , &
    !     feptype_1          , &
    !     zaerotype_bc1      , &
    !     zaerotype_bc2      , &
    !     zaerotype_dust1    , &
    !     zaerotype_dust2    , &
    !     zaerotype_dust3    , &
    !     zaerotype_dust4    , &
    !     ratio_C2N_diatoms  , &
    !     ratio_C2N_sp       , &
    !     ratio_C2N_phaeo    , &
    !     ratio_chl2N_diatoms, &
    !     ratio_chl2N_sp     , &
    !     ratio_chl2N_phaeo  , &
    !     F_abs_chl_diatoms  , &
    !     F_abs_chl_sp       , &
    !     F_abs_chl_phaeo    , &
    !     ratio_C2N_proteins

    use ice_colpkg, only: &
         colpkg_init_parameters

    type(domain_type), intent(inout) :: &
         domain

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
         config_ocean_heat_transfer_type, &
         config_sea_freezing_temperature_type, &
         config_skeletal_bgc_flux_type

    logical, pointer :: &
         config_calc_surface_temperature, &
         config_use_form_drag, &
         config_use_high_frequency_coupling, &
         config_use_ocean_mixed_layer, &
         config_calc_surface_stresses, &
         config_use_vertical_tracers, &
         config_scale_initial_vertical_bgc, &
         config_use_vertical_biochemistry, &
         config_use_shortwave_bioabsorption, &
         config_use_skeletal_biochemistry, &
         config_use_vertical_zsalinity, &
         config_use_modal_aerosols

    real(kind=RKIND), pointer :: &
         config_min_friction_velocity, &
         config_rapid_mode_channel_radius, &
         config_rapid_model_critical_Ra, &
         config_rapid_mode_aspect_ratio, &
         config_slow_mode_drainage_strength, &
         config_slow_mode_critical_porosity, &
         config_congelation_ice_porosity, &
         config_visible_ice_albedo, &
         config_infrared_ice_albedo, &
         config_visible_snow_albedo, &
         config_infrared_snow_albedo, &
         config_variable_albedo_thickness_limit, &
         config_ice_shortwave_tuning_parameter, &
         config_pond_shortwave_tuning_parameter, &
         config_snow_shortwave_tuning_parameter, &
         config_temp_change_snow_grain_radius_change, &
         config_max_melting_snow_grain_radius, &
         config_algae_absorption_coefficient, &
         config_ridiging_efolding_scale, &
         config_ratio_ridging_work_to_PE, &
         config_snow_to_ice_transition_depth, &
         config_pond_flushing_timescale, &
         config_min_meltwater_retained_fraction, &
         config_max_meltwater_retained_fraction, &
         config_pond_depth_to_fraction_ratio, &
         config_snow_on_pond_ice_tapering_parameter, &
         config_critical_pond_ice_thickness, &
         config_biogrid_bottom_molecular_sublayer, &
         config_bio_gravity_drainage_length_scale, &
         config_biogrid_top_molecular_sublayer, &
         config_new_ice_fraction_biotracer, &
         config_fraction_biotracer_in_frazil, &
         config_zsalinity_molecular_sublayer, &
         config_zsalinity_gravity_drainage_scale, &
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
         config_rapid_mobile_to_stationary_time, &
         config_long_mobile_to_stationary_time, &
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
         config_ratio_C_to_N_proteins

    integer, pointer :: &
         config_boundary_layer_iteration_number

    call MPAS_pool_get_config(domain % configs, "config_thermodynamics_type", config_thermodynamics_type)
    call MPAS_pool_get_config(domain % configs, "config_heat_conductivity_type", config_heat_conductivity_type)
    call MPAS_pool_get_config(domain % configs, "config_ocean_heat_transfer_type", config_ocean_heat_transfer_type)
    call MPAS_pool_get_config(domain % configs, "config_calc_surface_temperature", config_calc_surface_temperature)
    call MPAS_pool_get_config(domain % configs, "config_min_friction_velocity", config_min_friction_velocity)
    call MPAS_pool_get_config(domain % configs, "config_rapid_mode_channel_radius", config_rapid_mode_channel_radius)
    call MPAS_pool_get_config(domain % configs, "config_rapid_model_critical_Ra", config_rapid_model_critical_Ra)
    call MPAS_pool_get_config(domain % configs, "config_rapid_mode_aspect_ratio", config_rapid_mode_aspect_ratio)
    call MPAS_pool_get_config(domain % configs, "config_slow_mode_drainage_strength", config_slow_mode_drainage_strength)
    call MPAS_pool_get_config(domain % configs, "config_slow_mode_critical_porosity", config_slow_mode_critical_porosity)
    call MPAS_pool_get_config(domain % configs, "config_congelation_ice_porosity", config_congelation_ice_porosity)
    call MPAS_pool_get_config(domain % configs, "config_shortwave_type", config_shortwave_type)
    call MPAS_pool_get_config(domain % configs, "config_albedo_type", config_albedo_type)
    call MPAS_pool_get_config(domain % configs, "config_visible_ice_albedo", config_visible_ice_albedo)
    call MPAS_pool_get_config(domain % configs, "config_infrared_ice_albedo", config_infrared_ice_albedo)
    call MPAS_pool_get_config(domain % configs, "config_visible_snow_albedo", config_visible_snow_albedo)
    call MPAS_pool_get_config(domain % configs, "config_infrared_snow_albedo", config_infrared_snow_albedo)
    call MPAS_pool_get_config(domain % configs, "config_variable_albedo_thickness_limit", config_variable_albedo_thickness_limit)
    call MPAS_pool_get_config(domain % configs, "config_ice_shortwave_tuning_parameter", config_ice_shortwave_tuning_parameter)
    call MPAS_pool_get_config(domain % configs, "config_pond_shortwave_tuning_parameter", config_pond_shortwave_tuning_parameter)
    call MPAS_pool_get_config(domain % configs, "config_snow_shortwave_tuning_parameter", config_snow_shortwave_tuning_parameter)
    call MPAS_pool_get_config(domain % configs, "config_temp_change_snow_grain_radius_change", &
                                                 config_temp_change_snow_grain_radius_change)
    call MPAS_pool_get_config(domain % configs, "config_max_melting_snow_grain_radius", config_max_melting_snow_grain_radius)
    call MPAS_pool_get_config(domain % configs, "config_algae_absorption_coefficient", config_algae_absorption_coefficient)
    call MPAS_pool_get_config(domain % configs, "config_ice_strength_formulation", config_ice_strength_formulation)
    call MPAS_pool_get_config(domain % configs, "config_ridging_participation_function", config_ridging_participation_function)
    call MPAS_pool_get_config(domain % configs, "config_ridging_redistribution_function", config_ridging_redistribution_function)
    call MPAS_pool_get_config(domain % configs, "config_ridiging_efolding_scale", config_ridiging_efolding_scale)
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
    call MPAS_pool_get_config(domain % configs, "config_pond_flushing_timescale", config_pond_flushing_timescale)
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
    call MPAS_pool_get_config(domain % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)
    call MPAS_pool_get_config(domain % configs, "config_biogrid_bottom_molecular_sublayer", &
                                                 config_biogrid_bottom_molecular_sublayer)
    call MPAS_pool_get_config(domain % configs, "config_bio_gravity_drainage_length_scale", &
                                                 config_bio_gravity_drainage_length_scale)
    call MPAS_pool_get_config(domain % configs, "config_biogrid_top_molecular_sublayer", config_biogrid_top_molecular_sublayer)
    call MPAS_pool_get_config(domain % configs, "config_zsalinity_gravity_drainage_scale", config_zsalinity_gravity_drainage_scale)
    call MPAS_pool_get_config(domain % configs, "config_new_ice_fraction_biotracer", config_new_ice_fraction_biotracer)
    call MPAS_pool_get_config(domain % configs, "config_fraction_biotracer_in_frazil", config_fraction_biotracer_in_frazil)
    call MPAS_pool_get_config(domain % configs, "config_zsalinity_molecular_sublayer", config_zsalinity_molecular_sublayer)
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
    call MPAS_pool_get_config(domain % configs, "config_rapid_mobile_to_stationary_time", config_rapid_mobile_to_stationary_time)
    call MPAS_pool_get_config(domain % configs, "config_long_mobile_to_stationary_time", config_long_mobile_to_stationary_time)
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

    call colpkg_init_parameters(&
         config_cice_int("config_thermodynamics_type", config_thermodynamics_type), &
         config_heat_conductivity_type, &
         config_ocean_heat_transfer_type, &
         config_calc_surface_temperature, &
         config_min_friction_velocity, &
         config_rapid_mode_channel_radius, &
         config_rapid_model_critical_Ra, &
         config_rapid_mode_aspect_ratio, &
         config_slow_mode_drainage_strength, &
         config_slow_mode_critical_porosity, &
         config_congelation_ice_porosity, &
         config_shortwave_type, &
         config_albedo_type, &
         config_visible_ice_albedo, &
         config_infrared_ice_albedo, &
         config_visible_snow_albedo, &
         config_infrared_snow_albedo, &
         config_variable_albedo_thickness_limit, &
         config_ice_shortwave_tuning_parameter, &
         config_pond_shortwave_tuning_parameter, &
         config_snow_shortwave_tuning_parameter, &
         config_temp_change_snow_grain_radius_change, &
         config_max_melting_snow_grain_radius, &
         config_algae_absorption_coefficient, &
         config_cice_int("config_ice_strength_formulation", config_ice_strength_formulation), &
         config_cice_int("config_ridging_participation_function", config_ridging_participation_function), &
         config_cice_int("config_ridging_redistribution_function", config_ridging_redistribution_function), &
         config_ridiging_efolding_scale, &
         config_ratio_ridging_work_to_PE, &
         config_atmos_boundary_method, &
         config_calc_surface_stresses, &
         config_use_form_drag, &
         config_use_high_frequency_coupling, &
         config_boundary_layer_iteration_number, &
         config_use_ocean_mixed_layer, &
         config_sea_freezing_temperature_type, &
         config_cice_int("config_itd_conversion_type", config_itd_conversion_type), &
         config_cice_int("config_category_bounds_type", config_category_bounds_type), &
         config_snow_to_ice_transition_depth, &
         config_pond_refreezing_type, &
         config_pond_flushing_timescale, &
         config_min_meltwater_retained_fraction, &
         config_max_meltwater_retained_fraction, &
         config_pond_depth_to_fraction_ratio, &
         config_snow_on_pond_ice_tapering_parameter, &
         config_critical_pond_ice_thickness, &
         config_skeletal_bgc_flux_type, &
         config_use_vertical_tracers, &
         config_scale_initial_vertical_bgc, &
         config_use_vertical_biochemistry, &
         config_use_shortwave_bioabsorption, &
         config_use_modal_aerosols, &
         config_use_skeletal_biochemistry, &
         config_use_vertical_zsalinity, &
         config_biogrid_bottom_molecular_sublayer, &
         config_bio_gravity_drainage_length_scale, &
         config_biogrid_top_molecular_sublayer, &
         config_new_ice_fraction_biotracer, &
         config_fraction_biotracer_in_frazil, &
         config_zsalinity_molecular_sublayer, &
         config_zsalinity_gravity_drainage_scale, &
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
         config_rapid_mobile_to_stationary_time, &
         config_long_mobile_to_stationary_time, &
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
         config_ratio_C_to_N_proteins)

    !-----------------------------------------------------------------------
    ! Parameters for thermodynamics
    !-----------------------------------------------------------------------

    ! ktherm:
    ! type of thermodynamics
    ! 0 = 0-layer approximation
    ! 1 = Bitz and Lipscomb 1999
    ! 2 = mushy layer theory
    !ktherm = config_cice_int("config_thermodynamics_type", config_thermodynamics_type)

    ! conduct:
    ! 'MU71' or 'bubbly'
    !conduct = config_heat_conductivity_type

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

    ! phi_i_mushy:
    ! liquid fraction of congelation ice
    !phi_i_mushy = config_congelation_ice_porosity

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

    ! R_ice:
    ! sea ice tuning parameter; +1 > 1sig increase in albedo
    !R_ice = config_ice_shortwave_tuning_parameter

    ! R_pnd:
    ! ponded ice tuning parameter; +1 > 1sig increase in albedo
    !R_pnd = config_pond_shortwave_tuning_parameter

    ! R_snw:
    ! snow tuning parameter; +1 > ~.01 change in broadband albedo
    !R_snw = config_snow_shortwave_tuning_parameter

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
    !mu_rdg = config_ridiging_efolding_scale

    ! Cf
    ! ratio of ridging work to PE change in ridging (kstrength = 1)
    !Cf = config_ratio_ridging_work_to_PE

    !-----------------------------------------------------------------------
    ! Parameters for atmosphere
    !-----------------------------------------------------------------------

    ! atmbndy:
    ! atmo boundary method, 'default' ('ccsm3') or 'constant'
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
    !oceanmixed_ice = config_use_ocean_mixed_layer

    ! fbot_xfer_type:
    ! transfer coefficient type for ice-ocean heat flux
    !fbot_xfer_type = config_ocean_heat_transfer_type

    ! tfrz_option:
    ! form of ocean freezing temperature
    ! 'minus1p8' = -1.8 C
    ! 'linear_salt' = -depressT * sss
    ! 'mushy' conforms with ktherm=2
    !tfrz_option = config_sea_freezing_temperature_type

    !-----------------------------------------------------------------------
    ! Parameters for the ice thickness distribution
    !-----------------------------------------------------------------------

    ! kitd:
    ! type of itd conversions
    !   0 = delta function
    !   1 = linear remap
    !kitd = config_cice_int("config_itd_conversion_type", config_itd_conversion_type)

    ! kcatbound:
    !   0 = old category boundary formula
    !   1 = new formula giving round numbers
    !   2 = WMO standard
    !   3 = asymptotic formula
    !kcatbound = config_cice_int("config_category_bounds_type", config_category_bounds_type)

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
    ! alter e-folding time scale for flushing
    !dpscale = config_pond_flushing_timescale

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

    ! solve_zsal:
    ! if true, update salinity profile from solve_S_dt
    !solve_zsal = config_use_vertical_zsalinity

    ! modal_aero:
    ! if true, use modal aerosal optical properties
    ! only for use with tr_aero or tr_zaero
    !modal_aero = config_use_shortwave_bioabsorption

    ! grid_o:
    ! for bottom flux
    !grid_o = config_biogrid_bottom_molecular_sublayer

    ! l_sk:
    ! characteristic diffusive scale (zsalinity) (m)
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

    ! grid_oS:
    ! for bottom flux (zsalinity)
    !grid_oS = config_zsalinity_molecular_sublayer

    ! l_skS:
    ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)
    !l_skS = config_zsalinity_gravity_drainage_scale

  end subroutine init_column_package_configs

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

    COLUMN_ERROR_WRITE(trim(config_name) COMMA ' has invalid value')
    COLUMN_ERROR_WRITE(trim(config_name) COMMA ': ' COMMA trim(config_value))
    COLUMN_ERROR_WRITE('valid options: ' COMMA trim(valid_options))
    call MPAS_dmpar_global_abort("MPAS-seaice: "//trim(config_name)//' has invalid value')

  end subroutine config_error

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
!  init_column_non_activated_pointers
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 5th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_non_activated_pointers(domain)

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
         pkgTracerZSalinityActive


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
       call MPAS_pool_get_package(block % packages, "pkgTracerZSalinityActive", pkgTracerZSalinityActive)


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
       if (.not. pkgTracerZSalinityActive) then
          call set_stand_in_tracer_array(block, "verticalSalinity")
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

  end subroutine init_column_non_activated_pointers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  finalize_column_non_activated_pointers
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 29th October 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine finalize_column_non_activated_pointers(domain)

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
         pkgTracerZSalinityActive

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
       call MPAS_pool_get_package(block % packages, "pkgTracerZSalinityActive", pkgTracerZSalinityActive)

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
       if (.not. pkgTracerZSalinityActive) then
          call finalize_stand_in_tracer_array(block, "verticalSalinity")
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

  end subroutine finalize_column_non_activated_pointers

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
!> \author Adrian K. Turner, LANL
!> \date 3rd April 2015
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
!  cice_column_initial_air_drag_coefficient
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date
!> \details
!>
!
!-----------------------------------------------------------------------

  function cice_column_initial_air_drag_coefficient() result(airDragCoefficient)

    use ice_constants_colpkg, only: &
         vonkar, zref, iceruf

    real(kind=RKIND) :: airDragCoefficient

    airDragCoefficient = (vonkar/log(zref/iceruf)) &
                       * (vonkar/log(zref/iceruf)) ! atmo drag for RASM

  end function cice_column_initial_air_drag_coefficient

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_reinitialize_fluxes
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 31st August 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_column_reinitialize_fluxes(domain)

    type(domain_type) :: domain

    ! atmospheric fluxes
    call cice_column_reinitialize_atmospheric_fluxes(domain)

    ! oceanic fluxes
    call cice_column_reinitialize_oceanic_fluxes(domain)

  end subroutine cice_column_reinitialize_fluxes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_reinitialize_atmospheric_fluxes
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 31st August 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_column_reinitialize_atmospheric_fluxes(domain)

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
         config_use_column_package

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

       call MPAS_pool_get_config(block % configs, "config_use_column_package", config_use_column_package)

       if (config_use_column_package) then

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

       endif ! config_use_column_package

       block => block % next
    end do

  end subroutine cice_column_reinitialize_atmospheric_fluxes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_reinitialize_oceanic_fluxes
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 31st August 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_column_reinitialize_oceanic_fluxes(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         oceanFluxesPool

    real(kind=RKIND), dimension(:), pointer :: &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         oceanShortwaveFlux

    logical, pointer :: &
         config_use_column_package

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_column_package", config_use_column_package)

       if (config_use_column_package) then

          call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxesPool)

          call MPAS_pool_get_array(oceanFluxesPool, "oceanFreshWaterFlux", oceanFreshWaterFlux)
          call MPAS_pool_get_array(oceanFluxesPool, "oceanSaltFlux", oceanSaltFlux)
          call MPAS_pool_get_array(oceanFluxesPool, "oceanHeatFlux", oceanHeatFlux)
          call MPAS_pool_get_array(oceanFluxesPool, "oceanShortwaveFlux", oceanShortwaveFlux)

          oceanFreshWaterFlux(:) = 0.0_RKIND
          oceanSaltFlux(:)       = 0.0_RKIND
          oceanHeatFlux(:)       = 0.0_RKIND
          oceanShortwaveFlux(:)  = 0.0_RKIND

       endif ! config_use_column_package

       block => block % next
    end do

  end subroutine cice_column_reinitialize_oceanic_fluxes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_column_tracer_object_bio_tracer_number
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 14 September 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_column_tracer_object_for_biogeochemistry(domain, tracerObject)

    use ice_colpkg, only: colpkg_init_zbgc

    type(domain_type), intent(in) :: &
         domain

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    logical, pointer :: &
         config_use_brine, &
         config_use_vertical_zsalinity, &
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

    ! save tracer array size
    nTracers_temp = tracerObject % nTracers
    tracerObject % nTracers = tracerObject % nTracersNotBio

    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)
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
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxAerosolType", maxAerosolType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxAlgaeType", maxAlgaeType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxDOCType", maxDOCType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxDICType", maxDICType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxDONType", maxDONType)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "maxIronType", maxIronType)

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

    call colpkg_init_zbgc(&
         nBioLayers, &
         nIceLayers, &
         nSnowLayers, &
         nAlgae, &
         nzAerosols, &
         nDOC, &
         nDIC, &
         nDON, &
         nDissolvedIron, &
         nParticulateIron, &
         tracerObject % firstAncestorMask, &
         tracerObject % parentIndex, &
         tracerObject % ancestorNumber, &
         tracerObject % ancestorIndices, &
         tracerObject % nBioTracersShortwave, &
         config_use_brine, &
         tracerObject % index_brineFraction,&
         tracerObject % nTracers, &
         tracerObject % nBioTracers, &
         tracerObject % index_nitrateConc, &
         tracerObject % index_ammoniumConc, &
         tracerObject % index_silicateConc, &
         tracerObject % index_DMSConc, &
         tracerObject % index_nonreactiveConc, &
         tracerObject % index_verticalSalinity, &
         tracerObject % index_algaeConc, &
         tracerObject % index_algalCarbon, &
         tracerObject % index_algalChlorophyll, &
         tracerObject % index_DOCConc, &
         tracerObject % index_DONConc, &
         tracerObject % index_DICConc, &
         tracerObject % index_verticalAerosolsConc, &
         tracerObject % index_DMSPpConc, &
         tracerObject % index_DMSPdConc, &
         tracerObject % index_dissolvedIronConc, &
         tracerObject % index_particulateIronConc, &
         tracerObject % index_mobileFraction, &
         config_use_nitrate, &
         config_use_ammonium, &
         config_use_silicate, &
         config_use_DMS, &
         config_use_nonreactive, &
         config_use_vertical_zsalinity, &
         use_nitrogen, &
         config_use_carbon, &
         config_use_chlorophyll, &
         config_use_DON, &
         config_use_iron,&
         config_use_zaerosols, &
         tracerObject % index_verticalAerosolsConcShortwave, &
         tracerObject % index_chlorophyllShortwave, &
         tracerObject % index_algaeConcLayer, &
         tracerObject % index_nitrateConcLayer, &
         tracerObject % index_ammoniumConcLayer, &
         tracerObject % index_silicateConcLayer, &
         tracerObject % index_DMSConcLayer, &
         tracerObject % index_DMSPpConcLayer, &
         tracerObject % index_DMSPdConcLayer, &
         tracerObject % index_algalCarbonLayer, &
         tracerObject % index_algalChlorophyllLayer, &
         tracerObject % index_DICConcLayer, &
         tracerObject % index_DOCConcLayer, &
         tracerObject % index_nonreactiveConcLayer, &
         tracerObject % index_DONConcLayer, &
         tracerObject % index_dissolvedIronConcLayer, &
         tracerObject % index_particulateIronConcLayer, &
         tracerObject % index_verticalAerosolsConcLayer, &
         tracerObject % index_humicsConc, &
         tracerObject % index_humicsConcLayer, &
         config_use_humics, &
         config_use_vertical_zsalinity, &
         config_use_skeletal_biochemistry, &
         config_use_vertical_tracers, &
         config_use_shortwave_bioabsorption, &
         config_use_vertical_biochemistry, &
         config_fraction_biotracer_in_frazil, &
         config_new_ice_fraction_biotracer, &
         tracerObject % index_LayerIndexToDataArray, &
         tracerObject % index_LayerIndexToBioIndex, &
         tracerObject % nTracersNotBio, &
         maxAlgaeType, &
         maxDOCType, &
         maxDICType, &
         maxDONType, &
         maxIronType, &
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
         config_long_mobile_to_stationary_time)

    ! check calculated tracer array size
    if (nTracers_temp /= tracerObject % nTracers) then
       COLUMN_ERROR_WRITE("init_column_tracer_object Error: nTracers_temp, nTracers")
       WRITE (stderrUnit,*) nTracers_temp
       WRITE (stderrUnit,*) tracerObject % nTracers
       call MPAS_dmpar_global_abort("MPAS-seaice: tracer number error")
    endif

  end subroutine init_column_tracer_object_for_biogeochemistry

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

    use ice_colpkg, only: &
         colpkg_init_bgc, &
         colpkg_init_hbrine, &
         colpkg_init_zsalinity

    type(domain_type), intent(inout) :: domain

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         biogeochemistry, &
         ocean_coupling, &
         tracers

    logical, pointer :: &
         config_use_brine, &
         config_use_vertical_zsalinity, &
         config_use_vertical_tracers, &
         config_use_skeletal_biochemistry, &
         config_do_restart_zsalinity, &
         config_do_restart_bgc, &
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
         rayleighCriteriaReal

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
         nCellsSolve, &
         nIceLayers, &
         nBioLayers, &
         nBioLayersP1, &
         nBioLayersP2, &
         nCategories, &
         nShortwaveBio, &
         nZBGCTracers, &
         maxAerosolType, &
         maxAlgaeType, &
         maxDOCType, &
         maxDICType, &
         maxDONType, &
         maxIronType

    integer :: &
         iCell

    logical :: &
         abortFlag, &
         rayleighCriteria, &
         setGetPhysicsTracers, &
         setGetBGCTracers

    character(len=strKIND) :: &
         abortMessage

    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)
    call MPAS_pool_get_config(domain % configs, "config_do_restart_zsalinity", config_do_restart_zsalinity)
    call MPAS_pool_get_config(domain % configs, "config_do_restart_bgc", config_do_restart_bgc)
    call MPAS_pool_get_config(domain % configs, "config_do_restart_hbrine", config_do_restart_hbrine)
    call MPAS_pool_get_config(domain % configs, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)
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
       call MPAS_pool_get_array(biogeochemistry, "rayleighCriteriaReal", rayleighCriteriaReal)
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

       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
       call MPAS_pool_get_array(ocean_coupling, "seaSurfaceSalinity", seaSurfaceSalinity)

       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity, 1)
       call MPAS_pool_get_array(tracers, "brineFraction", brineFraction, 1)

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nCategories", nCategories)
       call MPAS_pool_get_dimension(block % dimensions, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(block % dimensions, "nBioLayers", nBioLayers)
       call MPAS_pool_get_dimension(block % dimensions, "nBioLayersP1", nBioLayersP1)
       call MPAS_pool_get_dimension(block % dimensions, "nBioLayersP2", nBioLayersP2)
       call MPAS_pool_get_dimension(block % dimensions, "nShortwaveBio", nShortwaveBio)
       call MPAS_pool_get_dimension(block % dimensions, "nZBGCTracers", nZBGCTracers)
       call MPAS_pool_get_dimension(block % dimensions, "maxAerosolType", maxAerosolType)
       call MPAS_pool_get_dimension(block % dimensions, "maxAlgaeType", maxAlgaeType)
       call MPAS_pool_get_dimension(block % dimensions, "maxDOCType", maxDOCType)
       call MPAS_pool_get_dimension(block % dimensions, "maxDICType", maxDICType)
       call MPAS_pool_get_dimension(block % dimensions, "maxDONType", maxDONType)
       call MPAS_pool_get_dimension(block % dimensions, "maxIronType", maxIronType)

       call colpkg_init_hbrine(&
            biologyGrid, &
            interfaceBiologyGrid, &
            verticalGrid, &
            interfaceGrid, &
            verticalShortwaveGrid, &
            nBioLayers, &
            nIceLayers, &
            config_snow_porosity_at_ice_surface)

       do iCell = 1, nCellsSolve

          if (.not. config_do_restart_hbrine) then
             ! initialize newly formed ice
             newlyFormedIce(:,iCell) = 1

             ! initialize brine fraction
             brineFraction(:,:,iCell) = 1.0_RKIND
          endif

          ! set the category tracer array
          call set_cice_tracer_array_category(block, tracerObject, iCell, setGetPhysicsTracers, setGetBGCTracers)

          if (config_use_vertical_zsalinity) then
             call colpkg_init_zsalinity(&
                  nBioLayers, &
                  tracerObject % nTracersNotBio, &
                  config_do_restart_zsalinity, &
                  rayleighCriteria, &
                  rayleighCriteriaReal(iCell), &
                  tracerObject % tracerArrayCategory(tracerObject % nTracersNotBio+1:tracerObject % nTracers,:), &
                  tracerObject % index_verticalSalinity, &
                  nCategories, &
                  seaSurfaceSalinity(iCell))
          endif

          if (config_use_vertical_tracers .or. config_use_skeletal_biochemistry) then
             call colpkg_init_bgc(&
                  config_dt, &
                  nCategories, &
                  nBioLayers, &
                  nIceLayers, &
                  tracerObject % nTracersNotBio, &
                  verticalGrid, &
                  interfaceBiologyGrid, &
                  config_do_restart_bgc, &
                  tracerObject % nTracers, &
                  tracerObject % nBiotracers, &
                  iceSalinity(:,:,iCell), &
                  tracerObject % tracerArrayCategory(tracerObject % nTracersNotBio+1:tracerObject % nTracers,:), &
                  seaSurfaceSalinity(iCell), &
                  oceanNitrateConc(iCell), &
                  oceanAmmoniumConc(iCell), &
                  oceanSilicateConc(iCell), &
                  oceanDMSPConc(iCell), &
                  oceanDMSConc(iCell), &
                  oceanAlgaeConc(:,iCell), &
                  oceanDOCConc(:,iCell), &
                  oceanDONConc(:,iCell), &
                  oceanDICConc(:,iCell), &
                  oceanDissolvedIronConc(:,iCell), &
                  oceanParticulateIronConc(:,iCell), &
                  oceanZAerosolConc(:,iCell), &
                  oceanHumicsConc(iCell), &
                  oceanBioConcentrations(:,iCell), &
                  maxAlgaeType, &
                  maxDOCType, &
                  maxDICType, &
                  maxDONType, &
                  maxIronType, &
                  nZBGCTracers, &
                  maxAerosolType, &
                  abortFlag, &
                  abortMessage)

             if (abortFlag) then
                 COLUMN_ERROR_WRITE("colpkg_init_bgc")
                 COLUMN_ERROR_WRITE(trim(abortMessage))
                 call MPAS_dmpar_global_abort("MPAS-seaice: "//trim(abortMessage))
             endif
          endif ! biogeochemistry

          ! get the category tracer array
          call get_cice_tracer_array_category(block, tracerObject, iCell, setGetPhysicsTracers, setGetBGCTracers)

       enddo ! iCell

       block => block % next
    end do

  end subroutine init_column_biogeochemistry_profiles

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_reinitialize_diagnostics_thermodynamics
!
!> \brief Reinitialize thermodynamics diagnostics
!> \author Adrian K. Turner, LANL
!> \date 27th September 2015
!> \details
!>  Reinitialize thermodynamics diagnostics
!
!-----------------------------------------------------------------------

  subroutine cice_column_reinitialize_diagnostics_thermodynamics(domain)

    use ice_constants_colpkg, only: &
         dragio

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         atmosFluxesPool, &
         meltGrowthRatesPool, &
         diagnosticsPool, &
         tracersAggregatePool, &
         pondsPool, &
         shortwavePool, &
         dragPool

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
         config_use_column_package

    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_package", config_use_column_package)

    if (config_use_column_package) then

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

          !fresh_ai  (:,:,:) = c0
          !fsalt_ai  (:,:,:) = c0
          !fhocn_ai  (:,:,:) = c0
          !fswthru_ai(:,:,:) = c0

          ! shortwave
          call MPAS_pool_get_subpool(block % structs, "shortwave", shortwavePool)

          call MPAS_pool_get_array(shortwavePool, "bareIceAlbedoCell", bareIceAlbedoCell)
          call MPAS_pool_get_array(shortwavePool, "snowAlbedoCell", snowAlbedoCell)
          call MPAS_pool_get_array(shortwavePool, "pondAlbedoCell", pondAlbedoCell)

          bareIceAlbedoCell = 0.0_RKIND
          snowAlbedoCell    = 0.0_RKIND
          pondAlbedoCell    = 0.0_RKIND

          ! form drag
          call MPAS_pool_get_subpool(block % structs, "drag", dragPool)

          call MPAS_pool_get_array(dragPool, "oceanDragCoefficient", oceanDragCoefficient)
          call MPAS_pool_get_array(dragPool, "airDragCoefficient", airDragCoefficient)

          oceanDragCoefficient = dragio
          airDragCoefficient   = cice_column_initial_air_drag_coefficient()

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

          block => block % next
       end do

    endif ! config_use_column_package

  end subroutine cice_column_reinitialize_diagnostics_thermodynamics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_reinitialize_diagnostics_dynamics
!
!> \brief Reinitialize dynamics diagnostics
!> \author Adrian K. Turner, LANL
!> \date 27th September 2015
!> \details
!>  Reinitialize dynamics diagnostics
!
!-----------------------------------------------------------------------

  subroutine cice_column_reinitialize_diagnostics_dynamics(domain)

    use ice_constants_colpkg, only: &
         dragio

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
         config_use_column_package, &
         config_use_velocity_solver

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_package", config_use_column_package)
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_velocity_solver", config_use_velocity_solver)
    call MPAS_pool_get_config(domain % blocklist % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)

    if (config_use_column_package) then

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

    endif ! config_use_column_package

  end subroutine cice_column_reinitialize_diagnostics_dynamics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_reinitialize_diagnostics_bgc
!
!> \brief Reinitialize BGC diagnostics
!> \author Adrian K. Turner, LANL
!> \date 27th September 2015
!> \details
!>  Reinitialize BGC diagnostics
!
!-----------------------------------------------------------------------

  subroutine cice_column_reinitialize_diagnostics_bgc(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         biogeochemistryPool

    ! biogeochemistry
    real(kind=RKIND), dimension(:), pointer :: &
         primaryProduction, &
         netSpecificAlgalGrowthRate, &
         netBrineHeight, &
         zSalinityFlux, &
         zSalinityGDFlux

    real(kind=RKIND), dimension(:,:), pointer :: &
         oceanBioFluxes, &
         atmosIceBioFluxes, &
         snowIceBioFluxes, &
         totalVerticalBiologyIce, &
         totalVerticalBiologySnow

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         bioTracerShortwave

    logical, pointer :: &
         config_use_column_biogeochemistry, &
         config_use_column_shortwave, &
         config_use_column_package

    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_package", config_use_column_package)

    if (config_use_column_package) then

       block => domain % blocklist
       do while (associated(block))

          ! biogeochemistry
          call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

          if (config_use_column_biogeochemistry) then

             call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistryPool)

             call MPAS_pool_get_array(biogeochemistryPool, "primaryProduction", primaryProduction)
             call MPAS_pool_get_array(biogeochemistryPool, "netSpecificAlgalGrowthRate", netSpecificAlgalGrowthRate)
             call MPAS_pool_get_array(biogeochemistryPool, "netBrineHeight", netBrineHeight)
             call MPAS_pool_get_array(biogeochemistryPool, "zSalinityFlux", zSalinityFlux)
             call MPAS_pool_get_array(biogeochemistryPool, "zSalinityGDFlux", zSalinityGDFlux)
             call MPAS_pool_get_array(biogeochemistryPool, "oceanBioFluxes", oceanBioFluxes)
             call MPAS_pool_get_array(biogeochemistryPool, "atmosIceBioFluxes", atmosIceBioFluxes)
             call MPAS_pool_get_array(biogeochemistryPool, "snowIceBioFluxes", snowIceBioFluxes)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalBiologyIce", totalVerticalBiologyIce)
             call MPAS_pool_get_array(biogeochemistryPool, "totalVerticalBiologySnow", totalVerticalBiologySnow)

             primaryProduction          = 0.0_RKIND
             netSpecificAlgalGrowthRate = 0.0_RKIND
             netBrineHeight             = 0.0_RKIND
             zSalinityFlux              = 0.0_RKIND
             zSalinityGDFlux            = 0.0_RKIND
             oceanBioFluxes             = 0.0_RKIND
             atmosIceBioFluxes          = 0.0_RKIND
             snowIceBioFluxes           = 0.0_RKIND
             totalVerticalBiologyIce    = 0.0_RKIND
             totalVerticalBiologySnow   = 0.0_RKIND

          endif

          call MPAS_pool_get_config(block % configs, "config_use_column_shortwave", config_use_column_shortwave)

          if (config_use_column_biogeochemistry .or. config_use_column_shortwave) then

             call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistryPool)
             call MPAS_pool_get_array(biogeochemistryPool, "bioTracerShortwave", bioTracerShortwave)
             bioTracerShortwave         = 0.0_RKIND

          endif

          block => block % next
       end do

    endif ! config_use_column_package

  end subroutine cice_column_reinitialize_diagnostics_bgc

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

!-----------------------------------------------------------------------
! Warning messages
!-----------------------------------------------------------------------

  subroutine column_write_warnings(warnings)

    character(len=strKINDWarnings), dimension(:), intent(in) :: &
         warnings

    integer :: &
         iWarning

    do iWarning = 1, size(warnings)

       write(stderrUnit,*) trim(warnings(iWarning))

    enddo ! iWarning

  end subroutine column_write_warnings

!-----------------------------------------------------------------------

end module cice_column
