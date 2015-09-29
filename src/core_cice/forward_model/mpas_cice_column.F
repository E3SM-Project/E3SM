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
       cice_column_reinitialize_fluxes

  ! constants
  real(kind=RKIND), parameter :: puny = 1e-11_RKIND
       
  ! tracer object
  type, private :: ciceTracerObjectType

     ! length of tracer array
     integer :: nTracers

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

     ! indexes of tracers in tracer array
     integer :: &
          index_surfaceTemperature, & ! nt_Tsfc
          index_iceEnthalpy, &        ! nt_qice
          index_snowEnthalpy, &       ! nt_qsno
          index_iceSalinity, &        ! nt_sice
          index_iceAge, &             ! nt_iage
          index_firstYearIceArea, &   ! nt_FY
          index_levelIceArea, &       ! nt_alvl
          index_levelIceVolume, &     ! nt_vlvl
          index_pondArea, &           ! nt_apnd
          index_pondDepth, &          ! nt_hpnd
          index_pondLidThickness, &   ! nt_ipnd
          index_aerosols              ! nt_aero

  end type ciceTracerObjectType

  type(ciceTracerObjectType), private :: ciceTracerObject

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

    ! set non activated variable pointers to other memory
    call init_column_non_activated_pointers(domain)

    ! initialize the column package tracer object
    call init_column_tracer_object(domain, ciceTracerObject)

    ! initialize the column package parameters
    call init_column_package_parameters(domain, ciceTracerObject)

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
         config_do_restart

    call MPAS_pool_get_config(domain % configs, "config_use_column_package", config_use_column_package)

    if (config_use_column_package) then

       ! initialize level ice tracers
       call init_column_level_ice_tracers(domain)
       
       ! initialize the itd thickness classes
       call init_column_itd(domain)
       
       ! initialize thermodynamic tracer profiles
       call init_column_thermodynamic_profiles(domain)
       
       ! history variables
       call init_column_history_variables(domain)

       ! shortwave
       call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)
       if (config_do_restart) &
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
    
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nCategories", nCategories)
    
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "initial", initial)

       call MPAS_pool_get_array(initial, "categoryThicknessLimits", categoryThicknessLimits)

       call colpkg_init_itd(&
            nCategories, &
            categoryThicknessLimits, &
            stdoutUnit)
    
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
             initialMeltingTemperatureProfile(iIceLayer,iCell) = colpkg_liquidus_temperature(initialSalinityProfileVertical(iIceLayer))

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

    use ice_colpkg, only: colpkg_init_orbit

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

    call MPAS_pool_get_config(domain % configs, "config_shortwave_type", config_shortwave_type)
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    if (trim(config_shortwave_type) == "dEdd") then

       abortFlag = .false.

       call colpkg_init_orbit(&
            stdoutUnit, &
            abortFlag, &
            abortMessage)

       if (abortFlag) then
          COLUMN_ERROR_WRITE("colpkg_init_orbit")
          COLUMN_ERROR_WRITE(trim(abortMessage))
          call MPAS_dmpar_global_abort(trim(abortMessage)) 
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

          shortwaveScalingFactor(iCell)   = 1.0_RKIND

          do iCategory = 1, nCategories

             ! aggregate albedos
             if (iceAreaCategory(1,iCategory,iCell) > puny) then

                albedoVisibleDirectCell(iCell)  = albedoVisibleDirectCell(iCell)  + &
                     albedoVisibleDirectCategory(iCategory,iCell)  * iceAreaCategory(1,iCategory,iCell)
                albedoVisibleDiffuseCell(iCell) = albedoVisibleDiffuseCell(iCell) + &
                     albedoVisibleDiffuseCategory(iCategory,iCell) * iceAreaCategory(1,iCategory,iCell)
                albedoIRDirectCell(iCell)       = albedoIRDirectCell(iCell)       + &
                     albedoIRDirectCategory(iCategory,iCell)       * iceAreaCategory(1,iCategory,iCell)
                albedoIRDiffuseCell(iCell)      = albedoIRDiffuseCell(iCell)      + &
                     albedoIRDiffuseCategory(iCategory,iCell)      * iceAreaCategory(1,iCategory,iCell)

                if (solarZenithAngleCosine(iCell) > puny) then ! sun above horizon

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
         config_use_level_ice

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         tracers

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         levelIceArea, &
         levelIceVolume

    call MPAS_pool_get_config(domain % configs, "config_use_level_ice", config_use_level_ice)

    if (config_use_level_ice) then
       
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
         !config_use_column_biogeochemistry, &
         config_use_column_itd_thermodynamics, &
         config_calc_surface_temperature

    real(kind=RKIND), pointer :: &
         config_dt

    call MPAS_pool_get_config(domain % configs, "config_use_column_package", config_use_column_package)
    
    if (config_use_column_package) then
    
       call MPAS_pool_get_config(domain % configs, "config_use_column_shortwave", config_use_column_shortwave)
       call MPAS_pool_get_config(domain % configs, "config_use_column_vertical_thermodynamics", config_use_column_vertical_thermodynamics)
       !call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
       call MPAS_pool_get_config(domain % configs, "config_use_column_itd_thermodynamics", config_use_column_itd_thermodynamics)
       call MPAS_pool_get_config(domain % configs, "config_calc_surface_temperature", config_calc_surface_temperature)

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

       !call mpas_timer_start("Column biogeochemistry")
       !if (config_use_column_biogeochemistry) &
       ! call column biogeochemistry ! not implemented yet
       !call mpas_timer_stop("Column biogeochemistry")

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

       call cice_debug_write_out_fields(domain, "cice_column_update_state", 1)
       call mpas_timer_start("Column update state")
       call cice_column_update_state(domain, config_dt, config_dt)
       call mpas_timer_stop("Column update state")
       call cice_debug_write_out_fields(domain, "cice_column_update_state", 2)

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
         config_use_column_ridging

    type(MPAS_pool_type), pointer :: &
         velocitySolver

    real(kind=RKIND), pointer :: &
         dynamicsTimeStep

    call MPAS_pool_get_config(domain % configs, "config_use_column_package", config_use_column_package)

    if (config_use_column_package) then
    
       call MPAS_pool_get_config(domain % configs, "config_use_column_ridging", config_use_column_ridging)

       call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolver)
       call MPAS_pool_get_array(velocitySolver, "dynamicsTimeStep", dynamicsTimeStep)

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

       call cice_debug_write_out_fields(domain, "cice_column_update_state", 3)
       call mpas_timer_start("Column update state")
       call cice_column_update_state(domain, dynamicsTimeStep, 0.0_RKIND)
       call mpas_timer_stop("Column update state")
       call cice_debug_write_out_fields(domain, "cice_column_update_state", 4)

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

    use ice_colpkg, only: colpkg_step_therm1

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

    call cice_debug_write_out_fields(domain, "column_vertical_thermodynamics", 1)

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

          abortFlag = .false.

          call colpkg_step_therm1(&
               config_dt, &
               nCategories, &
               nIceLayers, &
               nSnowLayers, &
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
               evaporativeHeatFlux(iCell), &
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
               abortMessage, &
               stdoutUnit)

          ! code abort
          if (abortFlag) then
             COLUMN_ERROR_WRITE("column_vertical_thermodynamics")
             COLUMN_ERROR_WRITE(trim(abortMessage))
             call MPAS_dmpar_global_abort(trim(abortMessage))
          endif

          ! aerosol
          if (config_use_aerosols) then

             do iCategory = 1, nCategories
                do iAerosol = 1, nAerosols

                   if (snowVolumeCategory(1,iCategory,iCell) > puny) &
                        specificSnowAerosol(iAerosol, :, iCategory) = specificSnowAerosol(iAerosol, :, iCategory) / snowVolumeCategory(1,iCategory,iCell)

                   if (iceVolumeCategory(1,iCategory,iCell) > puny) &
                        specificIceAerosol(iAerosol, :, iCategory)  = specificIceAerosol(iAerosol, :, iCategory)  / iceVolumeCategory(1,iCategory,iCell)

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

    call cice_debug_write_out_fields(domain, "column_vertical_thermodynamics", 2)
    
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

    use ice_colpkg, only: colpkg_step_therm2

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
         config_update_ocean_fluxes

    ! dimensions
    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nIceLayers, &
         nSnowLayers, &
         nAerosols

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
         categoryThicknessLimits

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

    ! local
    integer :: &
         iCell, &
         iCategory

    logical, dimension(:), allocatable :: &
         newlyFormedIceLogical

    logical :: &
         abortFlag

    character(len=strKIND) :: &
         abortMessage

    real(kind=RKIND) :: &
         dayOfYear

    call cice_debug_write_out_fields(domain, "column_itd_thermodynamics", 1)

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

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)
       call MPAS_pool_get_dimension(mesh, "nAerosols", nAerosols)

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

       call MPAS_pool_get_array(ponds, "pondFreshWaterFlux", pondFreshWaterFlux)

       call MPAS_pool_get_array(aerosols, "oceanAerosolFlux", oceanAerosolFlux)

       call MPAS_pool_get_array(biogeochemistry, "newlyFormedIce", newlyFormedIce)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioFluxes", oceanBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioConcentrations", oceanBioConcentrations)

       call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)
       call MPAS_pool_get_array(initial, "categoryThicknessLimits", categoryThicknessLimits)

       call MPAS_pool_get_array(diagnostics, "freezeOnset", freezeOnset)

       ! newly formed ice
       allocate(newlyFormedIceLogical(nCategories))
       
       do iCell = 1, nCellsSolve

          ! newly formed ice
          do iCategory = 1, nCategories
             newlyFormedIceLogical(iCategory) = (newlyFormedIce(iCategory,iCell) == 1)
          enddo ! iCategory

          ! set the category tracer array
          call set_cice_tracer_array_category(block, ciceTracerObject, iCell)

          abortFlag = .false.

          call colpkg_step_therm2(&
               config_dt, &
               nCategories, &
               nAerosols, &
               nIcelayers, &
               nSnowLayers, &
               categoryThicknessLimits, &
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
               oceanAerosolFlux(:,iCell), &
               newlyFormedIceLogical(:), &
               oceanBioFluxes(:,iCell), &
               oceanBioConcentrations(:,iCell), &
               abortFlag, &
               abortMessage, &
               stdoutUnit, &
               freezeOnset(iCell), &
               dayOfYear)

          ! get category tracer array
          call get_cice_tracer_array_category(block, ciceTracerObject, iCell)

          ! code abort
          if (abortFlag) then
             COLUMN_ERROR_WRITE("column_itd_thermodynamics")
             COLUMN_ERROR_WRITE(trim(abortMessage))
             call MPAS_dmpar_global_abort(trim(abortMessage))
          endif

       enddo ! iCell

       ! newly formed ice
       deallocate(newlyFormedIceLogical)

       block => block % next
    end do

    call cice_debug_write_out_fields(domain, "column_itd_thermodynamics", 2)

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

    call cice_debug_write_out_fields(domain, "column_prep_radiation", 1)

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

    call cice_debug_write_out_fields(domain, "column_prep_radiation", 2)

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

    use ice_colpkg, only: colpkg_step_radiation

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
         aerosols

    ! configs
    real(kind=RKIND), pointer :: &
         config_dt

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
         bareIceAlbedoCategory, &
         snowAlbedoCategory, &
         pondAlbedoCategory, &
         effectivePondAreaCategory, &
         pondSnowDepthDifference, &
         pondLidMeltFluxFraction, &
         aerosolMassExtinctionCrossSection, &
         aerosolSingleScatterAlbedo, &
         aerosolAsymmetryParameter, &
         albedoVisibleDirectCategory, &
         albedoVisibleDiffuseCategory, &
         albedoIRDirectCategory, &
         albedoIRDiffuseCategory

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
         iceBodyAerosol
    
    ! aerosols array
    real(kind=RKIND), dimension(:,:), allocatable :: &
         aerosolsArray

    ! local
    integer :: &
         iCell, &
         iCategory, &
         iAerosol

    real(kind=RKIND) :: &
         dayOfYear, &
         lonCellColumn

    integer :: &
         secondsIntoDay, &
         daysInYear

    ! CCSM coupled variables
    character(len=strKIND), parameter :: &
         CCSMCalendarType = ""

    real(kind=RKIND), parameter :: &
         CCSMDayNextShortwaveCalculation = 0.0_RKIND

    call cice_debug_write_out_fields(domain, "column_radiation", 1)

    ! day of year
    call get_day_of_year(clock, dayOfYear)

    ! seconds into day
    call get_seconds_into_day(clock, secondsIntoDay)

    ! get days in year
    call get_days_in_year(domain, clock, daysInYear)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "ponds", ponds)
       call MPAS_pool_get_subpool(block % structs, "aerosols", aerosols)

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)


       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)
       call MPAS_pool_get_dimension(mesh, "nAerosols", nAerosols)


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

       call MPAS_pool_get_array(atmos_coupling, "shortwaveVisibleDirectDown", shortwaveVisibleDirectDown)
       call MPAS_pool_get_array(atmos_coupling, "shortwaveVisibleDiffuseDown", shortwaveVisibleDiffuseDown)
       call MPAS_pool_get_array(atmos_coupling, "shortwaveIRDirectDown", shortwaveIRDirectDown)
       call MPAS_pool_get_array(atmos_coupling, "shortwaveIRDiffuseDown", shortwaveIRDiffuseDown)
       call MPAS_pool_get_array(atmos_coupling, "snowfallRate", snowfallRate)

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

       call MPAS_pool_get_array(ponds, "pondSnowDepthDifference", pondSnowDepthDifference)
       call MPAS_pool_get_array(ponds, "pondLidMeltFluxFraction", pondLidMeltFluxFraction)

       call MPAS_pool_get_array(aerosols, "aerosolMassExtinctionCrossSection", aerosolMassExtinctionCrossSection)
       call MPAS_pool_get_array(aerosols, "aerosolSingleScatterAlbedo", aerosolSingleScatterAlbedo)
       call MPAS_pool_get_array(aerosols, "aerosolAsymmetryParameter", aerosolAsymmetryParameter)

       ! aerosols array
       allocate(aerosolsArray(4*nAerosols,nCategories))

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
          if (lonCellColumn > pii) lonCellColumn = lonCellColumn - 2.0 * pii

          call colpkg_step_radiation(&
               config_dt, &
               nCategories, & 
               nIceLayers, &
               nSnowLayers, &
               nAerosols, &
               iceAreaCategory(1,:,iCell), &
               iceVolumeCategory(1,:,iCell), &
               snowVolumeCategory(1,:,iCell), & 
               surfaceTemperature(1,:,iCell), &
               levelIceArea(1,:,iCell), &
               pondArea(1,:,iCell), &
               pondDepth(1,:,iCell), &
               pondLidThickness(1,:,iCell), &
               aerosolsArray, &
               latCell(iCell), & 
               lonCellColumn, &
               CCSMCalendarType, &
               daysInYear, &
               CCSMDayNextShortwaveCalculation, &
               dayOfYear, &
               secondsIntoDay, &
               aerosolMassExtinctionCrossSection(:,:), &
               aerosolSingleScatterAlbedo(:,:), &
               aerosolAsymmetryParameter(:,:), &
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
               pondSnowDepthDifference(:,iCell), &
               pondLidMeltFluxFraction(:,iCell), &
               stdoutUnit, &
               .false., &
               lInitialization)

       enddo ! iCell

       ! aerosols array
       deallocate(aerosolsArray)

       block => block % next
    end do

    call cice_debug_write_out_fields(domain, "column_radiation", 2)

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

    use ice_colpkg, only: colpkg_step_ridge

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
         initial

    ! configs
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
         nAerosols

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

    ! local
    integer :: &
         iCell, &
         iCategory

    logical, dimension(:), allocatable :: &
         newlyFormedIceLogical

    logical :: &
         abortFlag

    character(len=strKIND) :: &
         abortMessage

    call cice_debug_write_out_fields(domain, "column_ridging", 1)

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

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_dynamics_subcycle_number", config_dynamics_subcycle_number)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nAerosols)

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

       call MPAS_pool_get_array(initial, "categoryThicknessLimits", categoryThicknessLimits)

       ! newly formed ice
       allocate(newlyFormedIceLogical(nCategories))

       do iCell = 1, nCellsSolve

          ! newly formed ice
          do iCategory = 1, nCategories
             newlyFormedIceLogical(iCategory) = (newlyFormedIce(iCategory,iCell) == 1)
          enddo ! iCategory

          ! set the category tracer array
          call set_cice_tracer_array_category(block, ciceTracerObject, iCell)

          abortFlag = .false.

          call colpkg_step_ridge(&
               config_dt, &        
               config_dynamics_subcycle_number, &
               nIceLayers, &
               nSnowLayers, &
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
               openingRateRidge(iCell),  &
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
               oceanBioFluxes(:,iCell), &
               stdoutUnit, &
               abortFlag, &  
               abortMessage)
          
          ! get category tracer array
          call get_cice_tracer_array_category(block, ciceTracerObject, iCell)

          ! code abort
          if (abortFlag) then
             COLUMN_ERROR_WRITE("column_ridging")
             COLUMN_ERROR_WRITE(trim(abortMessage))
             call MPAS_dmpar_global_abort(trim(abortMessage))
          endif

       enddo ! iCell

       ! newly formed ice
       deallocate(newlyFormedIceLogical)

       block => block % next
    end do

    call cice_debug_write_out_fields(domain, "column_ridging", 2)

  end subroutine column_ridging

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

  subroutine cice_column_update_state(domain, dt, iceAgeTimeOffset)

    type(domain_type), intent(inout) :: domain

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
         iceAreaTendencyTransport, &
         iceVolumeTendencyTransport, &
         iceAgeTendencyTransport

    integer, pointer :: &
         nCellsSolve
    
    integer :: &
         iCell

    logical, pointer :: &
         config_use_ice_age

    ! halo update of state variables
    call mpas_timer_start("Column halo update")
    call cice_column_state_halo_update(domain)
    call mpas_timer_stop("Column halo update")

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

       call MPAS_pool_get_array(diagnostics, "iceAreaTendencyTransport", iceAreaTendencyTransport)
       call MPAS_pool_get_array(diagnostics, "iceVolumeTendencyTransport", iceVolumeTendencyTransport)
       call MPAS_pool_get_array(diagnostics, "iceAgeTendencyTransport", iceAgeTendencyTransport)

       do iCell = 1, nCellsSolve

          iceAreaTendencyTransport(iCell) = &
               (iceAreaCell(iCell) - iceAreaTendencyTransport(iCell)) / dt
          iceVolumeTendencyTransport(iCell) = &
               (iceVolumeCell(iCell) - iceVolumeTendencyTransport(iCell)) / dt

          if (config_use_ice_age) then
             if (iceAgeTimeOffset > 0.0_RKIND) then
                
                if (iceAgeCell(iCell) > 0.0_RKIND) &
                     iceAgeTendencyTransport(iCell) = &
                          (iceAgeCell(iCell) - iceAgeTendencyTransport(iCell) - iceAgeTimeOffset) / dt

             else

                iceAgeTendencyTransport(iCell) = &
                     (iceAgeCell(iCell) - iceAgeTendencyTransport(iCell)) / dt
          
             endif
          endif

       enddo ! iCell

       block => block % next
    enddo

  end subroutine cice_column_update_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_state_halo_update
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 31st March 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_column_state_halo_update(domain)

    use mpas_dmpar

    type(domain_type), intent(inout) :: domain

    type(MPAS_pool_type), pointer :: &
         tracers

    type(MPAS_pool_iterator_type) :: &
         tracersIterator

    type(MPAS_pool_field_info_type) :: &
         tracerInfo
    
    type(field3DReal), pointer :: &
         tracerField3DReal

    ! get the tracers pool
    call MPAS_pool_get_subpool(domain % blocklist % structs, "tracers", tracers)

    ! iterate over the members of the tracers pool
    call MPAS_pool_begin_iteration(tracers)
    do while (MPAS_pool_get_next_member(tracers, tracersIterator))

       if (tracersIterator % memberType == MPAS_POOL_FIELD) then

          ! get the tracer info type
          call MPAS_pool_get_field_info(tracers, tracersIterator % memberName, tracerInfo)

          ! only exchange halos if tracer is active
          if (tracerInfo % isActive) then

             ! select tracer rank
             if (tracerInfo % nDims == 3) then

                ! get the field
                call MPAS_pool_get_field(tracers, tracersIterator % memberName, tracerField3DReal, 1)

                ! exchange the halo
                call MPAS_dmpar_exch_halo_field(tracerField3DReal)

             endif

          endif

       endif

    end do ! tracers pool iterator

  end subroutine cice_column_state_halo_update

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

          ! set the category tracer array
          call set_cice_tracer_array_category(block, ciceTracerObject, iCell)

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
          call get_cice_tracer_array_cell(block, ciceTracerObject, iCell)

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
         config_include_pond_freshwater_feedback

    type(MPAS_pool_type), pointer :: &
         shortwave, &
         atmosCoupling, &
         tracers, &
         ponds, &
         oceanFluxes

    real(kind=RKIND), dimension(:), pointer :: &
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
         oceanFreshWaterFlux

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

    real(kind=RKIND), pointer :: &
         config_dt

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    integer :: &
         iCell, &
         iCategory

    call cice_debug_write_out_fields(domain, "cice_column_coupling_prep", 1)

    call MPAS_pool_get_config(domain % configs, "config_use_ocean_mixed_layer", config_use_ocean_mixed_layer)
    call MPAS_pool_get_config(domain % configs, "config_dt", config_dt)
    call MPAS_pool_get_config(domain % configs, "config_include_pond_freshwater_feedback", config_include_pond_freshwater_feedback)

    if (config_use_ocean_mixed_layer) &
         call cice_column_ocean_mixed_layer(domain)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
       call MPAS_pool_get_subpool(block % structs, "ponds", ponds)
       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxes)

       call MPAS_pool_get_dimension(tracers, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(tracers, "nCategories", nCategories)

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

       do iCell = 1, nCellsSolve

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

          !-----------------------------------------------------------------
          ! Save net shortwave for scaling factor in shortwaveScalingFactor
          !-----------------------------------------------------------------

          shortwaveScalingFactor(iCell) = &
               shortwaveVisibleDirectDown(iCell)  * (1.0_RKIND - albedoVisibleDirectArea(iCell)) + &
               shortwaveVisibleDiffuseDown(iCell) * (1.0_RKIND - albedoVisibleDiffuseArea(iCell)) + &
               shortwaveIRDirectDown(iCell)       * (1.0_RKIND - albedoIRDirectArea(iCell)) + &
               shortwaveIRDiffuseDown(iCell)      * (1.0_RKIND - albedoIRDiffuseArea(iCell))
          
       enddo ! iCell


       block => block % next
    enddo

    call cice_column_scale_fluxes(domain)

    call cice_debug_write_out_fields(domain, "cice_column_coupling_prep", 2)

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
         oceanFluxes

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         airStressCellU, &
         airStressCellV, &
         sensibleHeatFlux, &
         latentHeatFlux, &
         absorbedShortwaveFlux, &
         longwaveUp, &
         evaporativeHeatFlux, &
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
         seaFreezingTemperature

    real(kind=RKIND) :: &
         iceAreaInverse

    integer, pointer :: &
         nCellsSolve

    integer :: &
         iCell

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregate)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolver)
       call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxes)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCoupling)
       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxes)

       call MPAS_pool_get_dimension(tracersAggregate, "nCellsSolve", nCellsSolve)

       call MPAS_pool_get_array(tracersAggregate, "iceAreaCell", iceAreaCell)

       call MPAS_pool_get_array(velocitySolver, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocitySolver, "airStressCellV", airStressCellV)

       call MPAS_pool_get_array(atmosFluxes, "sensibleHeatFlux", sensibleHeatFlux)
       call MPAS_pool_get_array(atmosFluxes, "latentHeatFlux", latentHeatFlux)
       call MPAS_pool_get_array(atmosFluxes, "evaporativeHeatFlux", evaporativeHeatFlux)
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

       do iCell = 1, nCellsSolve

          if (iceAreaCell(iCell) > 0.0_RKIND) then

             iceAreaInverse = 1.0_RKIND / iceAreaCell(iCell)

             airStressCellU(iCell)              = airStressCellU(iCell)              * iceAreaInverse
             airStressCellV(iCell)              = airStressCellV(iCell)              * iceAreaInverse
             sensibleHeatFlux(iCell)            = sensibleHeatFlux(iCell)            * iceAreaInverse
             latentHeatFlux(iCell)              = latentHeatFlux(iCell)              * iceAreaInverse
             absorbedShortwaveFlux(iCell)       = absorbedShortwaveFlux(iCell)       * iceAreaInverse
             longwaveUp(iCell)                  = longwaveUp(iCell)                  * iceAreaInverse
             evaporativeHeatFlux(iCell)         = evaporativeHeatFlux(iCell)         * iceAreaInverse
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

          else

             airStressCellU(iCell)              = 0.0_RKIND
             airStressCellV(iCell)              = 0.0_RKIND
             sensibleHeatFlux(iCell)            = 0.0_RKIND
             latentHeatFlux(iCell)              = 0.0_RKIND
             absorbedShortwaveFlux(iCell)       = 0.0_RKIND
             longwaveUp(iCell)                  = &
                  -ciceStefanBoltzmann * (seaFreezingTemperature(iCell) + ciceFreshWaterFreezingPoint)**4
             evaporativeHeatFlux(iCell)         = 0.0_RKIND
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
         airReferenceTemperature2mOcean, &
         airReferenceHumidity2mOcean, &
         longwaveUpOcean, &
         sensibleHeatFluxOcean, &
         latentHeatFluxOcean, &
         evaporativeHeatFluxOcean, &
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
       call MPAS_pool_get_array(oceanAtmosphere, "airReferenceTemperature2mOcean", airReferenceTemperature2mOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "airReferenceHumidity2mOcean", airReferenceHumidity2mOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "albedoVisibleDirectOcean", albedoVisibleDirectOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "albedoVisibleDiffuseOcean", albedoVisibleDiffuseOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "albedoIRDirectOcean", albedoIRDirectOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "albedoIRDiffuseOcean", albedoIRDiffuseOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "longwaveUpOcean", longwaveUpOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "sensibleHeatFluxOcean", sensibleHeatFluxOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "latentHeatFluxOcean", latentHeatFluxOcean)
       call MPAS_pool_get_array(oceanAtmosphere, "evaporativeHeatFluxOcean", evaporativeHeatFluxOcean)

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
               airReferenceTemperature2mOcean(iCell), &
               airReferenceHumidity2mOcean(iCell), &
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
               evaporativeHeatFluxOcean(iCell), &
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
         nCategories

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nCategories", nCategories)

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
         config_use_meltponds, &
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
    call MPAS_pool_get_config(domain % configs, "config_use_meltponds", config_use_meltponds)
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
    if (config_use_meltponds) then
       nTracers = nTracers + 1
       tracerObject % index_pondArea = nTracers
       nTracers = nTracers + 1
       tracerObject % index_pondDepth = nTracers
       if (config_use_level_meltponds) then
          nTracers = nTracers + 1
          tracerObject % index_pondLidThickness = nTracers
       endif
       if (config_use_topo_meltponds) then
          nTracers = nTracers + 1
          tracerObject % index_pondLidThickness = nTracers
       endif
    endif

    ! aerosols
    tracerObject % index_aerosols = indexMissingValue
    if (config_use_aerosols) then
       tracerObject % index_aerosols = nTracers + 1
    endif

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
         config_use_meltponds, &
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
    call MPAS_pool_get_config(domain % configs, "config_use_meltponds", config_use_meltponds)
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
       tracerObject % ancestorIndices(tracerObject % index_pondLidThickness,1) = tracerObject % index_levelIceArea  ! on level ice area
       
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
  
  subroutine set_cice_tracer_array_category(block, tracerObject, iCell)

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
         config_use_meltponds, &
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
    call MPAS_pool_get_config(block % configs, "config_use_meltponds", config_use_meltponds)
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
    if (config_use_meltponds) then
       tracerObject % tracerArrayCategory(nTracers,:) = pondArea(1,:,iCell)
       nTracers = nTracers + 1
       tracerObject % tracerArrayCategory(nTracers,:) = pondDepth(1,:,iCell)
       nTracers = nTracers + 1
    endif

    ! level or topo ponds
    if (config_use_meltponds .and. &
       (config_use_level_meltponds .or. config_use_topo_meltponds)) then
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

  subroutine get_cice_tracer_array_category(block, tracerObject, iCell)

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
         config_use_meltponds, &
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
    call MPAS_pool_get_config(block % configs, "config_use_meltponds", config_use_meltponds)
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
    if (config_use_meltponds) then
       pondArea(1,:,iCell) = tracerObject % tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
       pondDepth(1,:,iCell) = tracerObject % tracerArrayCategory(nTracers,:)
       nTracers = nTracers + 1
    endif

    ! level or topo ponds
    if (config_use_meltponds .and. &
       (config_use_level_meltponds .or. config_use_topo_meltponds)) then
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
  
  subroutine set_cice_tracer_array_cell(block, tracerObject, iCell)

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
         config_use_meltponds, &
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
    call MPAS_pool_get_config(block % configs, "config_use_meltponds", config_use_meltponds)
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
    if (config_use_meltponds) then
       tracerObject % tracerArrayCell(nTracers) = pondAreaCell(iCell)
       nTracers = nTracers + 1
       tracerObject % tracerArrayCell(nTracers) = pondDepthCell(iCell)
       nTracers = nTracers + 1
    endif

    ! level or topo ponds
    if (config_use_meltponds .and. &
       (config_use_level_meltponds .or. config_use_topo_meltponds)) then
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

  subroutine get_cice_tracer_array_cell(block, tracerObject, iCell)

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
         config_use_meltponds, &
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
    call MPAS_pool_get_config(block % configs, "config_use_meltponds", config_use_meltponds)
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
    if (config_use_meltponds) then
       pondAreaCell(iCell) = tracerObject % tracerArrayCell(nTracers)
       nTracers = nTracers + 1
       pondDepthCell(iCell) = tracerObject % tracerArrayCell(nTracers)
       nTracers = nTracers + 1
    endif

    ! level or topo ponds
    if (config_use_meltponds .and. &
       (config_use_level_meltponds .or. config_use_topo_meltponds)) then
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

  end subroutine get_cice_tracer_array_cell

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
         config_use_cesm_meltponds, &
         config_use_form_drag, &
         config_use_level_ice, &
         config_use_level_meltponds, &
         config_use_meltponds
    
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
    call MPAS_pool_get_config(domain % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_form_drag", config_use_form_drag)
    call MPAS_pool_get_config(domain % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_meltponds", config_use_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_ocean_heat_transfer_type", config_ocean_heat_transfer_type)
    call MPAS_pool_get_config(domain % configs, "config_sea_freezing_temperature_type", config_sea_freezing_temperature_type)

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
       call config_error("config_ridging_participation_function", config_ridging_participation_function, "'Thorndike75' or 'exponential'")
    endif

    ! check config_ridging_redistribution_function value
    if (.not. (trim(config_ridging_redistribution_function) == "Hibler80" .or. &
               trim(config_ridging_redistribution_function) == "exponential")) then
       call config_error("config_ridging_redistribution_function", config_ridging_redistribution_function, "'Hibler80' or 'exponential'")
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
       call config_error("config_category_bounds_type", config_category_bounds_type, "'single category', 'original', 'new', 'WMO' or 'asymptotic'")
    endif

    ! check config_pond_refreezing_type value 
    if (.not. (trim(config_pond_refreezing_type) == "cesm" .or. &
               trim(config_pond_refreezing_type) == "hlid")) then
       call config_error("config_pond_refreezing_type", config_pond_refreezing_type, "'cesm' or 'hlid'")
    endif

    !-----------------------------------------------------------------------
    ! Check combinations
    !-----------------------------------------------------------------------

    ! chech for itd remapping with only one category
    if (nCategories == 1 .and. trim(config_itd_conversion_type) == "linear remap") then
       COLUMN_ERROR_WRITE('Remapping the ITD is not allowed for nCategories=1')
       COLUMN_ERROR_WRITE("Use config_itd_conversion_type = 'delta function' with config_category_bounds_type = 'original'")
       COLUMN_ERROR_WRITE("or for column configurations use config_category_bounds_type = 'single category'")
       call MPAS_dmpar_global_abort('Incompatible nCategories and config_itd_conversion_type')
    endif

    ! check itd and category bounds discrepancy
    if (nCategories /= 1 .and. trim(config_category_bounds_type) == 'single category') then
       COLUMN_ERROR_WRITE("nCategories /= 1 .and. config_category_bounds_type = 'single category'")
       call MPAS_dmpar_global_abort('Incompatible nCategories and config_category_bounds_type')
    endif
    
    ! check config_snow_to_ice_transition_depth and level ponds
    if (config_use_level_meltponds .and. abs(config_snow_to_ice_transition_depth) > puny) then
       COLUMN_ERROR_WRITE("config_use_level_meltponds = .true. and config_snow_to_ice_transition_depth /= 0")
       call MPAS_dmpar_global_abort('Incompatible config_use_level_meltponds and config_snow_to_ice_transition_depth')
    endif

    ! check cesm ponds and freezing lids inconsistency
    if (config_use_cesm_meltponds .and. trim(config_pond_refreezing_type) /= "cesm") then
       COLUMN_ERROR_WRITE("config_use_cesm_meltponds = .true. and config_pond_refreezing_type /= 'cesm'")
       call MPAS_dmpar_global_abort('Incompatible config_use_cesm_meltponds and config_pond_refreezing_type')
    endif

    ! check dEdd shortwave if using ponds
    if (trim(config_shortwave_type) /= 'dEdd' .and. config_use_meltponds .and. config_calc_surface_temperature) then
       COLUMN_ERROR_WRITE("config_shortwave_type) /= 'dEdd' .and. config_use_meltponds = .true. .and. config_calc_surface_temperature ==.true.")
       call MPAS_dmpar_global_abort('Incompatible config_shortwave_type, config_use_meltponds and config_calc_surface_temperature')
    endif

    ! check range of config_min_meltwater_retained_fraction and config_max_meltwater_retained_fraction
    if (config_min_meltwater_retained_fraction < 0.0_RKIND .or. &
        config_min_meltwater_retained_fraction > 1.0_RKIND) then
       COLUMN_ERROR_WRITE('config_min_meltwater_retained_fraction out of bounds')
       call MPAS_dmpar_global_abort('config_min_meltwater_retained_fraction out of bounds')
    endif
    if (config_max_meltwater_retained_fraction < 0.0_RKIND .or. &
        config_max_meltwater_retained_fraction > 1.0_RKIND) then
       COLUMN_ERROR_WRITE('config_max_meltwater_retained_fraction out of bounds')
       call MPAS_dmpar_global_abort('config_max_meltwater_retained_fraction out of bounds')
    endif

    ! check not mushy physics and dont calculate surface temperature
    if (trim(config_thermodynamics_type) == "mushy" .and. .not. config_calc_surface_temperature) then
       COLUMN_ERROR_WRITE("config_thermodynamics_type = 'mushy' and config_calc_surface_temperature = .false.")
       call MPAS_dmpar_global_abort('Incompatible config_thermodynamics_type and config_calc_surface_temperature')
    endif

    ! check not form drag with constant atmosphere boundary method
    if (config_use_form_drag .and. trim(config_atmos_boundary_method) == "constant") then
       COLUMN_ERROR_WRITE("config_use_form_drag = .true. and config_atmos_boundary_method = 'constant'")
       call MPAS_dmpar_global_abort('Incompatible config_use_form_drag and config_atmos_boundary_method')
    endif

    ! check not form drag with not calculating surface stresses
    if (config_use_form_drag .and. .not. config_calc_surface_stresses) then
       COLUMN_ERROR_WRITE("config_use_form_drag = .true. and config_calc_surface_stresses = .false.")
       call MPAS_dmpar_global_abort('Incompatible config_use_form_drag and config_calc_surface_stresses')
    endif

    ! check am not using form drag with cesm ponds
    if (config_use_form_drag .and. config_use_cesm_meltponds) then
       COLUMN_ERROR_WRITE("config_use_form_drag = .true. and config_use_cesm_meltponds = .true.")
       call MPAS_dmpar_global_abort('Incompatible config_use_form_drag and config_use_cesm_meltponds')
    endif

    ! check using form drag but not level ice
    if (config_use_form_drag .and. .not. config_use_level_ice) then
       COLUMN_ERROR_WRITE("config_use_form_drag = .true. and config_use_level_ice = .false.")
       call MPAS_dmpar_global_abort('Incompatible config_use_form_drag and config_use_level_ice')
    endif

    ! check form drag and ocean heat flux type
    if (.not. config_use_form_drag .and. trim(config_ocean_heat_transfer_type) == "Cdn_ocn") then
       COLUMN_ERROR_WRITE("config_use_form_drag = .false. and config_ocean_heat_transfer_type == 'Cdn_ocn'")
       call MPAS_dmpar_global_abort('Incompatible config_use_form_drag and config_ocean_heat_transfer_type')
    endif

    ! check thermodynamic type and sea freezing temperature type
    if (trim(config_thermodynamics_type) == "BL99" .and. trim(config_sea_freezing_temperature_type) /= "linear_salt") then
       COLUMN_ERROR_WRITE("config_thermodynamics_type == 'BL99' and config_sea_freezing_temperature_type /= 'linear_salt'")
       call MPAS_dmpar_global_abort('Incompatible config_thermodynamics_type and config_sea_freezing_temperature_type')
    endif
    if (trim(config_thermodynamics_type) == "mushy" .and. trim(config_sea_freezing_temperature_type) /= "mushy") then
       COLUMN_ERROR_WRITE("config_thermodynamics_type == 'mushy' and config_sea_freezing_temperature_type /= 'mushy'")
       call MPAS_dmpar_global_abort('Incompatible config_thermodynamics_type and config_sea_freezing_temperature_type')
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
         config_use_meltponds, &
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols, &
         config_use_brine

    call MPAS_pool_get_config(domain % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(domain % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_meltponds", config_use_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_aerosols", config_use_aerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_brine", config_use_brine)

    call colpkg_init_tracer_flags(&
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_meltponds, &
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols, &
         config_use_brine)

    !tr_iage      = config_use_ice_age
    !tr_FY        = config_use_first_year_ice
    !tr_lvl       = config_use_level_ice
    !tr_pond      = config_use_meltponds
    !tr_pond_cesm = config_use_cesm_meltponds
    !tr_pond_lvl  = config_use_level_meltponds
    !tr_pond_topo = config_use_topo_meltponds
    !tr_aero      = config_use_aerosols
    !tr_brine     = config_use_brine

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
    !     nbtrcr
    use ice_colpkg, only: &
         colpkg_init_tracer_numbers

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    call colpkg_init_tracer_numbers(&
         tracerObject % nTracers, &
         0)

    !ntrcr = tracerObject % nTracers

    !nbtrcr = 0

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
    !     nt_Tsfc, & ! ice/snow temperature
    !     nt_qice, & ! volume-weighted ice enthalpy (in layers)
    !     nt_qsno, & ! volume-weighted snow enthalpy (in layers)
    !     nt_sice, & ! volume-weighted ice bulk salinity (CICE grid layers)
    !     nt_fbri, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
    !     nt_iage, & ! volume-weighted ice age
    !     nt_FY, & ! area-weighted first-year ice area
    !     nt_alvl, & ! level ice area fraction
    !     nt_vlvl, & ! level ice volume fraction
    !     nt_apnd, & ! melt pond area fraction
    !     nt_hpnd, & ! melt pond depth
    !     nt_ipnd, & ! melt pond refrozen lid thickness
    !     nt_aero, & ! starting index for aerosols in ice
    !     nt_bgc_N_sk, & ! algae (skeletal layer)
    !     nt_bgc_C_sk, & ! 
    !     nt_bgc_chl_sk, & ! 
    !     nt_bgc_Nit_sk, & ! nutrients (skeletal layer) 
    !     nt_bgc_Am_sk, & ! 
    !     nt_bgc_Sil_sk, & !
    !     nt_bgc_DMSPp_sk, & ! trace gases (skeletal layer)
    !     nt_bgc_DMSPd_sk, & ! 
    !     nt_bgc_DMS_sk, & ! 
    !     nt_bgc_Nit_ml, & ! nutrients (ocean mixed layer) 
    !     nt_bgc_Am_ml, & ! 
    !     nt_bgc_Sil_ml, & !
    !     nt_bgc_DMSP_ml, & ! trace gases (ocean mixed layer)
    !     nt_bgc_DMS_ml
    use ice_colpkg, only: &
         colpkg_init_tracer_indices

    type(ciceTracerObjectType), intent(in) :: &
         tracerObject

    call colpkg_init_tracer_indices(&
         tracerObject % index_surfaceTemperature, &
         tracerObject % index_iceEnthalpy, &
         tracerObject % index_snowEnthalpy, &
         tracerObject % index_iceSalinity, &
         0, &
         tracerObject % index_iceAge, &
         tracerObject % index_firstYearIceArea, &
         tracerObject % index_levelIceArea, &
         tracerObject % index_levelIceVolume, &
         tracerObject % index_pondArea, &
         tracerObject % index_pondDepth, &
         tracerObject % index_pondLidThickness, &
         tracerObject % index_aerosols, &
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    !nt_Tsfc = tracerObject % index_surfaceTemperature
    !nt_qice = tracerObject % index_iceEnthalpy
    !nt_qsno = tracerObject % index_snowEnthalpy
    !nt_sice = tracerObject % index_iceSalinity
    !nt_iage = tracerObject % index_iceAge
    !nt_FY   = tracerObject % index_firstYearIceArea
    !nt_alvl = tracerObject % index_levelIceArea
    !nt_vlvl = tracerObject % index_levelIceVolume
    !nt_apnd = tracerObject % index_pondArea
    !nt_hpnd = tracerObject % index_pondDepth
    !nt_ipnd = tracerObject % index_pondLidThickness
    !nt_aero = tracerObject % index_aerosols

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
         config_sea_freezing_temperature_type

    logical, pointer :: &
         config_calc_surface_temperature, &
         config_use_form_drag, &
         config_use_high_frequency_coupling, &
         config_use_ocean_mixed_layer, &
         config_calc_surface_stresses

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
         config_critical_pond_ice_thickness

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
    call MPAS_pool_get_config(domain % configs, "config_temp_change_snow_grain_radius_change", config_temp_change_snow_grain_radius_change)
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
    call MPAS_pool_get_config(domain % configs, "config_snow_on_pond_ice_tapering_parameter", config_snow_on_pond_ice_tapering_parameter)
    call MPAS_pool_get_config(domain % configs, "config_critical_pond_ice_thickness", config_critical_pond_ice_thickness)

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
         config_critical_pond_ice_thickness)
    
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
    ! critical Rayleigh number for rapid drainage mode
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
    call MPAS_dmpar_global_abort(trim(config_name)//' has invalid value')
    
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
         drag

    ! packages
    logical, pointer :: &
         pkgColumnTracerIceAgeActive, &
         pkgColumnTracerFirstYearIceActive, &
         pkgColumnTracerLevelIceActive, &
         pkgColumnTracerPondsActive, &
         pkgColumnTracerLidThicknessActive, &
         pkgColumnTracerAerosolsActive, &
         pkgColumnFormDragActive

    ! mesh stand-ins
    type(field1DReal), pointer :: &
         latCell ! nCells array

    ! drag variables
    type(field1DReal), pointer :: &
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

          call MPAS_pool_get_field(drag, "airOceanDragCoefficientRatio", airOceanDragCoefficientRatio)
          call MPAS_pool_get_field(drag, "oceanDragCoefficient", oceanDragCoefficient)
          call MPAS_pool_get_field(drag, "oceanDragCoefficientSkin", oceanDragCoefficientSkin)
          call MPAS_pool_get_field(drag, "oceanDragCoefficientFloe", oceanDragCoefficientFloe)
          call MPAS_pool_get_field(drag, "oceanDragCoefficientKeel", oceanDragCoefficientKeel)
          call MPAS_pool_get_field(drag, "airDragCoefficient", airDragCoefficient)
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

          airOceanDragCoefficientRatio % array => latCell % array
          oceanDragCoefficient % array         => latCell % array
          oceanDragCoefficientSkin % array     => latCell % array
          oceanDragCoefficientFloe % array     => latCell % array
          oceanDragCoefficientKeel % array     => latCell % array
          airDragCoefficient % array           => latCell % array
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
         evaporativeHeatFlux, &
         longwaveUp, &
         absorbedShortwaveFlux, &
         atmosReferenceTemperature2m, &
         atmosReferenceHumidity2m, &
         atmosReferenceSpeed10m
    
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxesPool)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwavePool)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCouplingPool)
       
       call MPAS_pool_get_array(velocitySolverPool, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressCellV", airStressCellV)

       call MPAS_pool_get_array(atmosFluxesPool, "sensibleHeatFlux", sensibleHeatFlux)
       call MPAS_pool_get_array(atmosFluxesPool, "latentHeatFlux", latentHeatFlux)
       call MPAS_pool_get_array(atmosFluxesPool, "evaporativeHeatFlux", evaporativeHeatFlux)
       call MPAS_pool_get_array(atmosFluxesPool, "longwaveUp", longwaveUp)

       call MPAS_pool_get_array(shortwavePool, "absorbedShortwaveFlux", absorbedShortwaveFlux)

       call MPAS_pool_get_array(atmosCouplingPool, "atmosReferenceTemperature2m", atmosReferenceTemperature2m)
       call MPAS_pool_get_array(atmosCouplingPool, "atmosReferenceHumidity2m", atmosReferenceHumidity2m)
       call MPAS_pool_get_array(atmosCouplingPool, "atmosReferenceSpeed10m", atmosReferenceSpeed10m)

       airStressCellU(:)              = 0.0_RKIND
       airStressCellV(:)              = 0.0_RKIND
       sensibleHeatFlux(:)            = 0.0_RKIND
       latentHeatFlux(:)              = 0.0_RKIND
       evaporativeHeatFlux(:)         = 0.0_RKIND
       longwaveUp(:)                  = 0.0_RKIND
       absorbedShortwaveFlux(:)       = 0.0_RKIND
       atmosReferenceTemperature2m(:) = 0.0_RKIND
       atmosReferenceHumidity2m(:)    = 0.0_RKIND
       atmosReferenceSpeed10m(:)      = 0.0_RKIND

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

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxesPool)

       call MPAS_pool_get_array(oceanFluxesPool, "oceanFreshWaterFlux", oceanFreshWaterFlux)
       call MPAS_pool_get_array(oceanFluxesPool, "oceanSaltFlux", oceanSaltFlux)
       call MPAS_pool_get_array(oceanFluxesPool, "oceanHeatFlux", oceanHeatFlux)
       call MPAS_pool_get_array(oceanFluxesPool, "oceanShortwaveFlux", oceanShortwaveFlux)

       oceanFreshWaterFlux(:) = 0.0_RKIND
       oceanSaltFlux(:)       = 0.0_RKIND
       oceanHeatFlux(:)       = 0.0_RKIND
       oceanShortwaveFlux(:)  = 0.0_RKIND

       block => block % next
    end do

  end subroutine cice_column_reinitialize_oceanic_fluxes

!-----------------------------------------------------------------------

end module cice_column

