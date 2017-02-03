!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_initialize
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module cice_initialize

  use mpas_dmpar
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_stream_manager
  use cice_debug

  implicit none

  private
  save

  public :: &
       cice_init, &
       cice_init_post_clock_advance

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_init(&
       domain, &
       clock, &
       dt)!{{{

    use cice_mesh, only: cice_init_mesh
    use cice_advection, only: cice_init_advection
    use cice_velocity_solver, only: cice_init_velocity_solver
    use cice_column, only: &
         cice_init_column_physics_package_parameters, &
         cice_init_column_physics_package_variables
    use cice_forcing, only: cice_forcing_init, cice_reset_coupler_fluxes

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    real(kind=RKIND), intent(in) :: &
         dt !< Input:

    type(block_type), pointer :: &
         block

    ! initialize junk values
    call init_junk_values(domain)

    ! initialize landice mask if needed for testing
    call init_test_ice_shelf_mask(domain)

    ! mesh initialization
    call cice_init_mesh(domain)

    ! init the basic column physics package
    call cice_init_column_physics_package_parameters(domain)

    ! init coupler fluxes
    call initialize_coupler_fields(domain)

    ! initialize forcing
    call cice_forcing_init(domain)

    ! init dynamics
    call cice_init_velocity_solver(domain)

    ! init advection
    call cice_init_advection(domain)

    ! column physics initialization
    call cice_init_column_physics_package_variables(domain, clock)

    ! init ice state
    call init_ice_state(domain, clock)

    ! initial halo exchange
    call initial_halo_exchanges(domain)

    ! coupler flux reset
    call cice_reset_coupler_fluxes(domain)

  end subroutine cice_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_post_clock_advance
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_init_post_clock_advance(&
       domain, &
       clock)!{{{

    use cice_column, only: &
         cice_init_column_shortwave

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_shortwave, &
         config_do_restart

    call MPAS_pool_get_config(domain % configs, "config_use_column_package", config_use_column_package)
    call MPAS_pool_get_config(domain % configs, "config_use_column_shortwave", config_use_column_shortwave)
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    if (config_use_column_package .and. config_use_column_shortwave .and. .not. config_do_restart) &
         call cice_init_column_shortwave(domain, clock)

  end subroutine cice_init_post_clock_advance

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_junk_values
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_junk_values(domain)!{{{

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh

    integer :: &
         njunk

    real(kind=RKIND), parameter :: &
         value_junk = -1.0e34_RKIND

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

       call MPAS_pool_get_dimension(mesh, "nCells", nCells)

       call MPAS_pool_get_array(mesh, "areaCell", areaCell)

       njunk = nCells + 1

       areaCell(njunk) = value_junk

       block => block % next
      end do

  end subroutine init_junk_values!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_state
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_state(&
       domain, &
       clock)!{{{

    use cice_testing, only: &
         cice_init_square_test_case_hex, &
         cice_init_square_point_test_case_hex
    use cice_column, only: &
         cice_column_aggregate

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_clock_type), intent(in) :: &
         clock

    type(block_type), pointer :: &
         block

    character(len=strKIND), pointer :: &
         config_initial_condition_type, &
         config_initial_velocity_type

    integer, dimension(:), pointer :: &
         interiorVertex

    type(MPAS_pool_type), pointer :: &
         icestate, &
         tracers_aggregate

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCellInitial, &
         iceAreaCell

    logical, pointer :: &
         config_do_restart, &
         config_use_column_package

    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    if (.not. config_do_restart) then

       call MPAS_pool_get_config(domain % configs, "config_initial_condition_type", config_initial_condition_type)
       call MPAS_pool_get_config(domain % configs, "config_initial_velocity_type", config_initial_velocity_type)

       ! set volumes/areas
       block => domain % blocklist
       do while (associated(block))

          if (trim(config_initial_condition_type) == "uniform") then

             call init_ice_state_uniform_ice(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "circle") then

             call init_ice_state_circle_of_ice(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "uniform_interior") then

             call init_ice_state_uniform_interior(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "special") then

             call init_ice_state_special(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "random_coverage") then

             call init_ice_state_random_coverage(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "square") then

             call cice_init_square_test_case_hex(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "square_point") then

             call cice_init_square_point_test_case_hex(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "cice_default") then

             call init_ice_cice_default(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "no_ice") then

             ! no ice initial condition
             continue

          else if (trim(config_initial_condition_type) == "none") then

             ! for test cases
             continue

          else

             write(stderrUnit,*) "Ice initial condition type unknown!"
             write(stderrUnit,*) "config_initial_condition_type: ", trim(config_initial_condition_type)
             call MPAS_dmpar_global_abort('MPAS-seaice: Ice initial condition type unknown')

          endif

          block => block % next
       end do

       ! halo exchanges
       call MPAS_dmpar_field_halo_exch(domain, 'iceAreaCategory', timeLevel=1)
       call MPAS_dmpar_field_halo_exch(domain, 'iceVolumeCategory', timeLevel=1)
       call MPAS_dmpar_field_halo_exch(domain, 'snowVolumeCategory', timeLevel=1)
       call MPAS_dmpar_field_halo_exch(domain, 'surfaceTemperature', timeLevel=1)

       ! other initialize variables
       block => domain % blocklist
       do while (associated(block))
          call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
          call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
          call MPAS_pool_get_array(icestate, "iceAreaCellInitial", iceAreaCellInitial)
          call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)
          iceAreaCellInitial = iceAreaCell
          block => block % next
       enddo

       ! set velocities
       block => domain % blocklist
       do while (associated(block))

          if (trim(config_initial_velocity_type) == "uniform") then

             call init_ice_velocity_uniform(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_velocity_type) == "ocean") then

             call init_ice_velocity_ocean(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_velocity_type) == "random") then

             call init_ice_velocity_random(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_velocity_type) == "special") then

             call init_ice_velocity_special(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_velocity_type) == "none") then

             continue

          else

             write(stderrUnit,*) "Ice initial velocity type unknown!"
             write(stderrUnit,*) "config_initial_velocity_type: ", trim(config_initial_velocity_type)
             call MPAS_dmpar_global_abort('MPAS-seaice: Ice initial velocity type unknown')

          endif

          block => block % next
       end do

       ! halo exchanges
       call MPAS_dmpar_field_halo_exch(domain, 'uVelocity')
       call MPAS_dmpar_field_halo_exch(domain, 'vVelocity')

    endif

    ! aggregate tracers even if restart
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_package", config_use_column_package)
    if (config_use_column_package) &
         call cice_column_aggregate(domain)

  end subroutine init_ice_state!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_state_uniform_ice
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_state_uniform_ice(&
       block, &
       configs)!{{{

    use cice_constants, only: &
         ciceDegreesToRadians

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate

    integer, pointer :: &
         nCellsSolve

    real(kind=RKIND), dimension(:), pointer :: &
         latCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         surfaceTemperatureCell

    real(kind=RKIND), pointer :: &
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume, &
         config_initial_latitude_north, &
         config_initial_latitude_south

    integer :: &
         iCell

    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)

    call MPAS_pool_get_dimension(mesh, "nCells", nCellsSolve)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)
    call MPAS_pool_get_array(tracers_aggregate, "surfaceTemperatureCell", surfaceTemperatureCell)

    do iCell = 1, nCellsSolve

       if (latCell(iCell) > config_initial_latitude_north * ciceDegreesToRadians .or. &
           latCell(iCell) < config_initial_latitude_south * ciceDegreesToRadians) then

          ! has ice
          iceAreaCategory(1,:,iCell)    = config_initial_ice_area
          iceVolumeCategory(1,:,iCell)  = config_initial_ice_volume
          snowVolumeCategory(1,:,iCell) = config_initial_snow_volume
          surfaceTemperature(1,:,iCell) = -1.0_RKIND

       else

          ! no ice
          iceAreaCategory(1,:,iCell)    = 0.0_RKIND
          iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
          snowVolumeCategory(1,:,iCell) = 0.0_RKIND
          surfaceTemperature(1,:,iCell) = 0.0_RKIND

       endif

    enddo ! iCell

    ! integrated quantities
    do iCell = 1, nCellsSolve

       iceAreaCell(iCell) = sum(iceAreaCategory(1,:,iCell))
       surfaceTemperatureCell(iCell) = -20.15_RKIND

    enddo ! iCell

  end subroutine init_ice_state_uniform_ice!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_cice_default
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date March 2nd 1015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_cice_default(&
       block, &
       configs)!{{{

    use cice_constants, only: &
         ciceDegreesToRadians

    use ice_colpkg, only: &
         colpkg_init_trcr, &
         colpkg_enthalpy_snow

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         atmos_coupling, &
         ocean_coupling, &
         initial

    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nIceLayers, &
         nSnowLayers

    real(kind=RKIND), dimension(:), pointer :: &
         latCell, &
         airTemperature, &
         seaSurfaceTemperature, &
         seaFreezingTemperature

      integer, dimension(:), pointer :: &
           landIceMask

    real(kind=RKIND), dimension(:,:), pointer :: &
         initialSalinityProfile, &
         initialMeltingTemperatureProfile

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         surfaceTemperature, &
         iceEnthalpy, &
         snowEnthalpy, &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         iceSalinity

    real(kind=RKIND), pointer :: &
         config_initial_latitude_north, &
         config_initial_latitude_south

    real(kind=RKIND), dimension(:), allocatable :: &
         initialCategoryIceArea, &
         initialCategoryIceThickness

    integer :: &
         iCell, &
         iCategory, &
         iIceLayer

    real(kind=RKIND), parameter :: &
         initialCategorySnowThickness = 0.2_RKIND

    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
    call MPAS_pool_get_subpool(block % structs, "initial", initial)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
    call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
    call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity, 1)

    call MPAS_pool_get_array(atmos_coupling, "airTemperature", airTemperature)

    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTemperature", seaSurfaceTemperature)
    call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)
    call MPAS_pool_get_array(ocean_coupling, "landIceMask", landIceMask)

    call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)
    call MPAS_pool_get_array(initial, "initialMeltingTemperatureProfile", initialMeltingTemperatureProfile)

    ! initial volumes and areas
    allocate(initialCategoryIceArea(nCategories))
    allocate(initialCategoryIceThickness(nCategories))

    call initial_category_areas_and_volumes(&
         block, &
         initialCategoryIceArea, &
         initialCategoryIceThickness)

    do iCell = 1, nCellsSolve

       if (seaSurfaceTemperature(iCell) <= seaFreezingTemperature(iCell) + 0.2_RKIND .and. &
          (latCell(iCell) > config_initial_latitude_north * ciceDegreesToRadians .or. &
           latCell(iCell) < config_initial_latitude_south * ciceDegreesToRadians) .and. &
          landIceMask(iCell) == 0) then

          ! has ice
          do iCategory = 1, nCategories

             iceAreaCategory(1,iCategory,iCell)    = initialCategoryIceArea(iCategory)
             iceVolumeCategory(1,iCategory,iCell)  = initialCategoryIceArea(iCategory) * initialCategoryIceThickness(iCategory)
             snowVolumeCategory(1,iCategory,iCell) = min(iceAreaCategory(1,iCategory,iCell) * initialCategorySnowThickness, &
                                                         0.2_RKIND * iceVolumeCategory(1,iCategory,iCell))

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

       else

          ! no ice
          iceAreaCategory(1,:,iCell)    = 0.0_RKIND
          iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
          snowVolumeCategory(1,:,iCell) = 0.0_RKIND

          surfaceTemperature(1,:,iCell) = seaFreezingTemperature(iCell)

          snowEnthalpy(1,:,iCell) = colpkg_enthalpy_snow(0.0_RKIND)

       endif

    enddo ! iCell

    ! clean up
    deallocate(initialCategoryIceArea)
    deallocate(initialCategoryIceThickness)

    ! get vertical salinity profiles
    call MPAS_pool_get_subpool(block % structs, "initial", initial)
    call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)

    do iCell = 1, nCellsSolve
       do iCategory = 1, nCategories
          do iIceLayer = 1, nIceLayers
             iceSalinity(iIceLayer,iCategory,iCell) = initialSalinityProfile(iIceLayer,iCell)
          enddo ! iIceLayer
       enddo ! iCategory
    enddo ! iCell

  end subroutine init_ice_cice_default!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  initial_category_areas_and_volumes
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date March 2nd 1015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine initial_category_areas_and_volumes(&
       block, &
       initialCategoryIceArea, &
       initialCategoryIceThickness)

    use cice_constants, only: &
         cicePuny

    use ice_colpkg, only: &
         colpkg_init_itd

    ! Note: the resulting average ice thickness
    ! tends to be less than hbar due to the
    ! nonlinear distribution of ice thicknesses

    type(block_type), intent(in) :: block

    real(kind=RKIND), dimension(:), intent(out) :: &
         initialCategoryIceArea, &
         initialCategoryIceThickness

    type(MPAS_pool_type), pointer :: &
         mesh, &
         initial

    real(kind=RKIND), dimension(:), pointer :: &
         categoryThicknessLimits

    real(kind=RKIND) :: &
         areaCategorySum

    integer, pointer :: &
         nCategories

    integer :: &
         iCategory

    logical, pointer :: &
         config_use_column_package

    real(kind=RKIND), parameter :: &
         thicknessWithLargestArea = 3.0_RKIND ! initial ice thickness with greatest area

    call MPAS_pool_get_subpool(block % structs, "initial", initial)
    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

    call MPAS_pool_get_array(initial, "categoryThicknessLimits", categoryThicknessLimits)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    !--------------------------------------------------------
    ! init caregory limits if needed
    !--------------------------------------------------------

    call MPAS_pool_get_config(block % configs, "config_use_column_package", config_use_column_package)

    if (.not. config_use_column_package) then

       call colpkg_init_itd(&
            nCategories, &
            categoryThicknessLimits, &
            stdoutUnit)

    endif

    !--------------------------------------------------------
    ! volumes
    !--------------------------------------------------------

    ! loop over categories to set initial volumes
    do iCategory = 1, nCategories-1

       ! middle of the category limits
       initialCategoryIceThickness(iCategory) = &
            0.5_RKIND * (categoryThicknessLimits(iCategory) + categoryThicknessLimits(iCategory+1))

    enddo ! iCategory

    ! thicknest category
    initialCategoryIceThickness(nCategories) = &
         categoryThicknessLimits(nCategories) + 1.0_RKIND

    !--------------------------------------------------------
    ! areas
    !--------------------------------------------------------

    ! initialize sum of areas in categories
    areaCategorySum = 0.0_RKIND

    do iCategory = 1, nCategories

       ! parabola, max at h=thicknessWithLargestArea, zero at h=0, 2*thicknessWithLargestArea
       initialCategoryIceArea(iCategory) = &
            max(0.0_RKIND, &
            (2.0_RKIND * thicknessWithLargestArea * initialCategoryIceThickness(iCategory) - &
             initialCategoryIceThickness(iCategory)**2))

       areaCategorySum = areaCategorySum + initialCategoryIceArea(iCategory)

    enddo ! iCategory

    ! normalize
    do iCategory = 1, nCategories

       initialCategoryIceArea(iCategory) = initialCategoryIceArea(iCategory) / &
            (areaCategorySum + cicePuny / nCategories)

    enddo ! iCategory

  end subroutine initial_category_areas_and_volumes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_state_uniform_interior_ice
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_state_uniform_interior(&
       block, &
       configs)!{{{

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         boundary

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         latCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    integer, dimension(:), pointer :: &
         interiorCell

    real(kind=RKIND), pointer :: &
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume, &
         config_initial_latitude_north, &
         config_initial_latitude_south

    integer :: &
         iCell

    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

    call MPAS_pool_get_array(boundary, "interiorCell", interiorCell)

    do iCell = 1, nCells

       if (interiorCell(iCell) == 1) then

          ! has ice
          iceAreaCategory(1,:,iCell)    = config_initial_ice_area
          iceVolumeCategory(1,:,iCell)  = config_initial_ice_volume
          snowVolumeCategory(1,:,iCell) = config_initial_snow_volume

       else

          ! no ice
          iceAreaCategory(1,:,iCell)    = 0.0_RKIND
          iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
          snowVolumeCategory(1,:,iCell) = 0.0_RKIND

       endif

    enddo ! iCell

  end subroutine init_ice_state_uniform_interior!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_state_circle_of_ice
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_state_circle_of_ice(&
       block, &
       configs)!{{{

    use cice_constants, only: &
         earthRadius

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         zCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature

    real(kind=RKIND), pointer :: &
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume

    logical, pointer :: &
         config_rotate_cartesian_grid

    integer :: &
         iCell

    real(kind=RKIND), parameter :: &
         circle_radius = 0.1_RKIND*earthRadius

    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, 'config_rotate_cartesian_grid', config_rotate_cartesian_grid)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    if (config_rotate_cartesian_grid) then

       do iCell = 1, nCells

          ! Put a circle of ice at the North Pole (which will be rotated to the grid equator)
          if (zCell(iCell) > 0.0_RKIND .and. sqrt(xCell(iCell)**2 + yCell(iCell)**2) < circle_radius) then

             ! has ice
             iceAreaCategory(1,:,iCell)    = config_initial_ice_area
             iceVolumeCategory(1,:,iCell)  = config_initial_ice_volume
             snowVolumeCategory(1,:,iCell) = config_initial_snow_volume
             surfaceTemperature(1,:,iCell) = 1.0_RKIND

          else

             ! no ice
             iceAreaCategory(1,:,iCell)    = 0.0_RKIND
             iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
             snowVolumeCategory(1,:,iCell) = 0.0_RKIND
             surfaceTemperature(1,:,iCell) = 0.0_RKIND

          endif

       enddo ! iCell

    else   ! no rotation

       do iCell = 1, nCells

          ! Put a circle of ice on the equator
          if (xCell(iCell) > 0.0_RKIND .and. sqrt(yCell(iCell)**2 + zCell(iCell)**2) < circle_radius) then   ! Greenwich meridian
!          if (xCell(iCell) < 0.0_RKIND .and. sqrt(yCell(iCell)**2 + zCell(iCell)**2) < circle_radius) then  ! Date Line

             ! has ice
             iceAreaCategory(1,:,iCell)    = config_initial_ice_area
             iceVolumeCategory(1,:,iCell)  = config_initial_ice_volume
             snowVolumeCategory(1,:,iCell) = config_initial_snow_volume
             surfaceTemperature(1,:,iCell) = 1.0_RKIND

          else

             ! no ice
             iceAreaCategory(1,:,iCell)    = 0.0_RKIND
             iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
             snowVolumeCategory(1,:,iCell) = 0.0_RKIND
             surfaceTemperature(1,:,iCell) = 0.0_RKIND

          endif

       enddo ! iCell

    endif ! config_rotate_cartesian_grid

  end subroutine init_ice_state_circle_of_ice!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_state_uniform_ice
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 30th July 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_state_special(&
       block, &
       configs)!{{{

    use cice_constants, only: &
         ciceRadiansToDegrees

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         latCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    real(kind=RKIND), pointer :: &
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume, &
         config_initial_latitude_north, &
         config_initial_latitude_south

    integer :: &
         iCell1, &
         iCell2, &
         iCell0, &
         iEdge1, &
         iEdge2, &
         iEdgeOnCell0, &
         iEdgeOnCell1, &
         iEdgeOnCell2

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         edgesOnCell

    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

    iceAreaCategory = 0.0_RKIND
    iceVolumeCategory = 0.0_RKIND
    snowVolumeCategory = 0.0_RKIND

    !iceAreaCategory = 1.0e-20_RKIND
    !iceVolumeCategory = 1.0_RKIND
    !snowVolumeCategory = 1.0_RKIND


    iCell1 = 2050
    iCell0 = 2051
    iCell2 = 2052

    do iEdgeOnCell0 = 1, nEdgesOnCell(iCell0)
       do iEdgeOnCell1 = 1, nEdgesOnCell(iCell1)
          if (edgesOnCell(iEdgeOnCell0,iCell0) == edgesOnCell(iEdgeOnCell1,iCell1)) then
             iEdge1 = edgesOnCell(iEdgeOnCell0,iCell0)
          endif
       enddo ! iEdgeOnCell1
    enddo ! iEdgeOnCell0

    do iEdgeOnCell0 = 1, nEdgesOnCell(iCell0)
       do iEdgeOnCell2 = 1, nEdgesOnCell(iCell2)
          if (edgesOnCell(iEdgeOnCell0,iCell0) == edgesOnCell(iEdgeOnCell2,iCell2)) then
             iEdge2 = edgesOnCell(iEdgeOnCell0,iCell0)
          endif
       enddo ! iEdgeOnCell1
    enddo ! iEdgeOnCell0

    iceAreaCategory(:,:,iCell1) = 1.0_RKIND
    iceVolumeCategory(:,:,iCell1) = 1.0_RKIND

    iceAreaCategory(:,:,iCell2) = 1.0_RKIND
    iceVolumeCategory(:,:,iCell2) = 1.0_RKIND

  end subroutine init_ice_state_special!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_state_uniform_ice
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 30th July 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_state_random_coverage(&
       block, &
       configs)!{{{

    use cice_constants, only: &
         ciceRadiansToDegrees

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         latCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature

    real(kind=RKIND), pointer :: &
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume, &
         config_initial_latitude_north, &
         config_initial_latitude_south

    integer :: &
         iCell

    real(kind=RKIND) :: &
         random, &
         random2

    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    do iCell = 1, nCells

       call random_number(random)

       if (random > 0.5_RKIND) then

          call random_number(random2)

          ! has ice
          iceAreaCategory(1,:,iCell)    = 10.0_RKIND**(-11.0_RKIND * random2)!config_initial_ice_area
          iceVolumeCategory(1,:,iCell)  = 10.0_RKIND**(-11.0_RKIND * random2)!config_initial_ice_volume
          snowVolumeCategory(1,:,iCell) = 10.0_RKIND**(-11.0_RKIND * random2)!config_initial_snow_volume
          surfaceTemperature(1,:,iCell) = 1.0_RKIND

       else

          ! no ice
          iceAreaCategory(1,:,iCell)    = 0.0_RKIND
          iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
          snowVolumeCategory(1,:,iCell) = 0.0_RKIND
          surfaceTemperature(1,:,iCell) = 0.0_RKIND

       endif

    enddo ! iCell

  end subroutine init_ice_state_random_coverage!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_velocity_uniform
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_velocity_uniform(&
       block, &
       configs)!{{{

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         velocity_solver, &
         boundary, &
         ocean_coupling

    integer, pointer :: &
         nVertices

    integer, dimension(:), pointer :: &
         interiorVertex, &
         landIceMaskVertex

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), pointer :: &
         config_initial_uvelocity, &
         config_initial_vvelocity

    integer :: &
         iVertex

    call MPAS_pool_get_config(configs, "config_initial_uvelocity", config_initial_uvelocity)
    call MPAS_pool_get_config(configs, "config_initial_vvelocity", config_initial_vvelocity)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)

    call MPAS_pool_get_array(ocean_coupling, "landIceMaskVertex", landIceMaskVertex)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    do iVertex = 1, nVertices

       if (interiorVertex(iVertex) == 1 .and. &
           landIceMaskVertex(iVertex) == 0) then

          uVelocity(iVertex) = config_initial_uvelocity
          vVelocity(iVertex) = config_initial_vvelocity

       else

          uVelocity(iVertex) = 0.0_RKIND
          vVelocity(iVertex) = 0.0_RKIND

       endif

    enddo ! iVertex

  end subroutine init_ice_velocity_uniform!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_velocity_ocean
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_velocity_ocean(&
       block, &
       configs)!{{{

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         velocity_solver, &
         boundary, &
         ocean_coupling

    integer, pointer :: &
         nVertices

    integer, dimension(:), pointer :: &
         interiorVertex, &
         landIceMaskVertex

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         uOceanVelocity, &
         vOceanVelocity

    integer :: &
         iVertex

    real(kind=RKIND), dimension(:), allocatable :: &
         uOceanVelocityVertex, &
         vOceanVelocityVertex

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)

    call MPAS_pool_get_array(ocean_coupling, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(ocean_coupling, "vOceanVelocity", vOceanVelocity)
    call MPAS_pool_get_array(ocean_coupling, "landIceMaskVertex", landIceMaskVertex)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    allocate(uOceanVelocityVertex(nVertices))
    allocate(vOceanVelocityVertex(nVertices))

    call cice_interpolate_cell_to_vertex(mesh, &
         uOceanVelocityVertex, &
         uOceanVelocity)

    call cice_interpolate_cell_to_vertex(mesh, &
         vOceanVelocityVertex, &
         vOceanVelocity)

    do iVertex = 1, nVertices

       if (interiorVertex(iVertex) == 1 .and. &
           landIceMaskVertex(iVertex) == 0) then

          uVelocity(iVertex) = uOceanVelocityVertex(iVertex)
          vVelocity(iVertex) = vOceanVelocityVertex(iVertex)

       else

          uVelocity(iVertex) = 0.0_RKIND
          vVelocity(iVertex) = 0.0_RKIND

       endif

    enddo ! iVertex

    deallocate(uOceanVelocityVertex)
    deallocate(vOceanVelocityVertex)

  end subroutine init_ice_velocity_ocean!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_velocity_random
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_velocity_random(&
       block, &
       configs)!{{{

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         velocity_solver, &
         boundary

    integer, pointer :: &
         nVertices

    integer, dimension(:), pointer :: &
         interiorVertex

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), pointer :: &
         config_initial_uvelocity, &
         config_initial_vvelocity

    integer :: &
         iVertex

    call MPAS_pool_get_config(configs, "config_initial_uvelocity", config_initial_uvelocity)
    call MPAS_pool_get_config(configs, "config_initial_vvelocity", config_initial_vvelocity)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call random_seed()

    do iVertex = 1, nVertices

       if (interiorVertex(iVertex) == 1) then

          call random_number(uVelocity(iVertex))
          call random_number(vVelocity(iVertex))

          uVelocity(iVertex) = config_initial_uvelocity * (uVelocity(iVertex) * 2.0_RKIND - 1.0_RKIND)
          vVelocity(iVertex) = config_initial_vvelocity * (vVelocity(iVertex) * 2.0_RKIND - 1.0_RKIND)

       else

          uVelocity(iVertex) = 0.0_RKIND
          vVelocity(iVertex) = 0.0_RKIND

       endif

    enddo ! iVertex

  end subroutine init_ice_velocity_random!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_velocity_special
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_velocity_special(&
       block, &
       configs)!{{{

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         velocity_solver, &
         boundary

    integer, pointer :: &
         nVertices

    integer, dimension(:), pointer :: &
         interiorVertex

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), pointer :: &
         config_initial_uvelocity, &
         config_initial_vvelocity

    integer :: &
         iVertex

    integer :: &
         iCell1, &
         iCell2, &
         iCell0, &
         iEdge1, &
         iEdge2, &
         iEdgeOnCell0, &
         iEdgeOnCell1, &
         iEdgeOnCell2

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         edgesOnCell

    call MPAS_pool_get_config(configs, "config_initial_uvelocity", config_initial_uvelocity)
    call MPAS_pool_get_config(configs, "config_initial_vvelocity", config_initial_vvelocity)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)

    uVelocity = 0.0_RKIND
    vVelocity = 0.0_RKIND

    iCell1 = 2050
    iCell0 = 2051
    iCell2 = 2052

    do iEdgeOnCell0 = 1, nEdgesOnCell(iCell0)
       do iEdgeOnCell1 = 1, nEdgesOnCell(iCell1)
          write(*,*) iEdgeOnCell0, iEdgeOnCell1, edgesOnCell(iEdgeOnCell0,iCell0), edgesOnCell(iEdgeOnCell1,iCell1)
          if (edgesOnCell(iEdgeOnCell0,iCell0) == edgesOnCell(iEdgeOnCell1,iCell1)) then
             iEdge1 = edgesOnCell(iEdgeOnCell0,iCell0)
          endif
       enddo ! iEdgeOnCell1
    enddo ! iEdgeOnCell0

    do iEdgeOnCell0 = 1, nEdgesOnCell(iCell0)
       do iEdgeOnCell2 = 1, nEdgesOnCell(iCell2)
          if (edgesOnCell(iEdgeOnCell0,iCell0) == edgesOnCell(iEdgeOnCell2,iCell2)) then
             iEdge2 = edgesOnCell(iEdgeOnCell0,iCell0)
          endif
       enddo ! iEdgeOnCell1
    enddo ! iEdgeOnCell0

    uVelocity(:) = 1.0_RKIND

  end subroutine init_ice_velocity_special!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  initial_halo_exchanges
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date Januray 13th 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine initial_halo_exchanges(domain)

    type(domain_type), intent(inout) :: &
         domain

    logical, pointer :: &
         pkgWeakActive

    call MPAS_pool_get_package(domain % blocklist % packages, "pkgWeakActive", pkgWeakActive)

    if (pkgWeakActive) then

       call MPAS_dmpar_field_halo_exch(domain, 'normalVectorTriangle')

    endif

  end subroutine initial_halo_exchanges

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
! initialize_coupler_fields
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 14th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine initialize_coupler_fields(domain)

    use ice_colpkg, only: &
         colpkg_liquidus_temperature, &
         colpkg_init_ocean_conc

    use cice_constants, only: &
         ciceStefanBoltzmann, &
         ciceFreshWaterFreezingPoint

    use cice_column, only: &
         cice_column_initial_air_drag_coefficient

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         atmosCoupling, &
         atmosForcing, &
         oceanCoupling, &
         atmosFluxes, &
         drag, &
         biogeochemistry, &
         aerosols

    real(kind=RKIND), dimension(:), pointer :: &
         uAirVelocity, &
         vAirVelocity, &
         seaSurfaceSalinity, &
         seaFreezingTemperature, &
         seaSurfaceTemperature, &
         windSpeed, &
         airDragCoefficient, &
         longwaveUp

    real(kind=RKIND), dimension(:), pointer :: &
         oceanNitrateConc, &
         oceanSilicateConc, &
         oceanAmmoniumConc, &
         oceanDMSConc, &
         oceanDMSPConc, &
         oceanHumicsConc, &
         carbonToNitrogenRatioAlgae, &
         carbonToNitrogenRatioDON

    real(kind=RKIND), dimension(:,:), pointer :: &
         atmosBioFluxes, &
         atmosBlackCarbonFlux, &
         atmosDustFlux, &
         oceanAlgaeConc, &
         oceanDOCConc, &
         oceanDICConc, &
         oceanDONConc, &
         oceanParticulateIronConc, &
         oceanDissolvedIronConc, &
         oceanZAerosolConc, &
         atmosAerosolFlux

    integer, pointer :: &
         maxAlgaeType, &
         maxDOCType, &
         maxDICType, &
         maxDONType, &
         maxIronType, &
         maxBCType, &
         maxDustType, &
         maxAerosolType, &
         nZBGCTracers

    logical, pointer :: &
         config_do_restart, &
         config_use_column_biogeochemistry, &
         config_use_column_package, &
         config_use_aerosols

    integer, pointer :: &
         nCells

    integer :: &
         iCell, &
         iBio, &
         indexj

    block => domain % blocklist
    do while (associated(block))

       !-------------------------------------------------------------
       ! Physics fluxes received from atmosphere
       !-------------------------------------------------------------

       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
       call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmosForcing)

       call MPAS_pool_get_array(atmosCoupling, "uAirVelocity", uAirVelocity)
       call MPAS_pool_get_array(atmosCoupling, "vAirVelocity", vAirVelocity)

       call MPAS_pool_get_array(atmosForcing, "windSpeed", windSpeed)

       windSpeed = sqrt(uAirVelocity**2 + vAirVelocity**2)

       !-------------------------------------------------------------
       ! Physics fluxes received from ocean
       !-------------------------------------------------------------

       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCoupling)
       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

       call MPAS_pool_get_dimension(mesh, "nCells", nCells)

       call MPAS_pool_get_array(oceanCoupling, "seaSurfaceSalinity", seaSurfaceSalinity)
       call MPAS_pool_get_array(oceanCoupling, "seaFreezingTemperature", seaFreezingTemperature)
       call MPAS_pool_get_array(oceanCoupling, "seaSurfaceTemperature", seaSurfaceTemperature)

       call MPAS_pool_get_config(block % configs, "config_do_restart", config_do_restart)

       do iCell = 1, nCells
          seaFreezingTemperature(iCell) = colpkg_liquidus_temperature(seaSurfaceSalinity(iCell))
       enddo ! iCell

       ! sea surface temperature is not initialized if we're restarting
       if (.not. config_do_restart) then
          seaSurfaceTemperature = seaFreezingTemperature
       endif

       !-------------------------------------------------------------
       ! fluxes sent to atmosphere
       !-------------------------------------------------------------

       call MPAS_pool_get_config(block % configs, "config_use_column_package", config_use_column_package)

       if (config_use_column_package) then

          call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxes)
          call MPAS_pool_get_array(atmosFluxes, "longwaveUp", longwaveUp)
          longwaveUp = -ciceStefanBoltzmann * ciceFreshWaterFreezingPoint**4

          call MPAS_pool_get_subpool(block % structs, "drag", drag)
          call MPAS_pool_get_array(drag, "airDragCoefficient", airDragCoefficient)
          airDragCoefficient = cice_column_initial_air_drag_coefficient()

       endif

       !-------------------------------------------------------------
       ! Bio fluxes received from atmosphere
       !-------------------------------------------------------------

       call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)

       if (config_use_aerosols) then
          call MPAS_pool_get_subpool(block % structs, "aerosols", aerosols)
          call MPAS_pool_get_array(aerosols, "atmosAerosolFlux", atmosAerosolFlux)

          atmosAerosolFlux = 1.e-12_RKIND
       endif

       call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

       if (config_use_column_biogeochemistry) then

          call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)

          call MPAS_pool_get_array(biogeochemistry, "atmosBioFluxes", atmosBioFluxes)
          call MPAS_pool_get_array(biogeochemistry, "atmosBlackCarbonFlux", atmosBlackCarbonFlux)
          call MPAS_pool_get_array(biogeochemistry, "atmosDustFlux", atmosDustFlux)

          call MPAS_pool_get_dimension(block % dimensions, "maxAlgaeType", maxAlgaeType)
          call MPAS_pool_get_dimension(block % dimensions, "maxDOCType", maxDOCType)
          call MPAS_pool_get_dimension(block % dimensions, "maxDICType", maxDICType)
          call MPAS_pool_get_dimension(block % dimensions, "maxDONType", maxDONType)
          call MPAS_pool_get_dimension(block % dimensions, "maxIronType", maxIronType)
          call MPAS_pool_get_dimension(block % dimensions, "maxBCType", maxBCType)
          call MPAS_pool_get_dimension(block % dimensions, "maxDustType", maxDustType)
          call MPAS_pool_get_dimension(block % dimensions, "maxAerosolType", maxAerosolType)
          call MPAS_pool_get_dimension(block % dimensions, "nZBGCTracers", nZBGCTracers)

          atmosBioFluxes = 0.0_RKIND
          do iCell = 1, nCells
             do iBio = 1, maxBCType
                atmosBlackCarbonFlux(iBio,iCell) =  1.e-12_RKIND
             enddo
             do iBio = 1, maxDustType
                atmosDustFlux(iBio,iCell) =  1.e-13_RKIND
             enddo
          enddo ! iCell

          !-------------------------------------------------------------
          ! Bio Concentrations received from ocean
          !-------------------------------------------------------------

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
          call MPAS_pool_get_array(biogeochemistry, "carbonToNitrogenRatioAlgae", carbonToNitrogenRatioAlgae)
          call MPAS_pool_get_array(biogeochemistry, "carbonToNitrogenRatioDON", carbonToNitrogenRatioDON)

          do iCell = 1, nCells

             call colpkg_init_ocean_conc(&
                  oceanAmmoniumConc(iCell), &
                  oceanDMSPConc(iCell), &
                  oceanDMSConc(iCell), &
                  oceanAlgaeConc(:,iCell), &
                  oceanDOCConc(:,iCell), &
                  oceanDICConc(:,iCell), &
                  oceanDONConc(:,iCell), &
                  oceanDissolvedIronConc(:,iCell), &
                  oceanParticulateIronConc(:,iCell), &
                  oceanHumicsConc(iCell), &
                  oceanNitrateConc(iCell), &
                  oceanSilicateConc(iCell),&
                  oceanZAerosolConc(:,iCell), &
                  maxDICType, &
                  maxDONType, &
                  maxIronType, &
                  maxAerosolType, &
                  carbonToNitrogenRatioAlgae, &
                  carbonToNitrogenRatioDON)

          enddo ! iCell

       endif ! config_use_column_biogeochemistry

       block => block % next
    enddo

  end subroutine initialize_coupler_fields

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
! init_test_ice_shelf_mask
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 9th April 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_test_ice_shelf_mask(domain)

    use cice_constants, only: &
         ciceDegreesToRadians, &
         ciceRadiansToDegrees

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         oceanCouplingPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         landIceMask, &
         latCell, &
         lonCell

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    real(kind=RKIND) :: &
         lonLimitWest = -66.10_RKIND, &
         lonLimitEast = -11.26_RKIND , &
         latLimit     = -71.94_RKIND

    logical, pointer :: &
         config_use_test_ice_shelf

    call MPAS_pool_get_config(domain % configs, "config_use_test_ice_shelf", config_use_test_ice_shelf)

    if (config_use_test_ice_shelf) then

       lonLimitWest = lonLimitWest + 360.0_RKIND
       lonLimitEast = lonLimitEast + 360.0_RKIND

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)
          call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

          call MPAS_pool_get_array(oceanCouplingPool, "landIceMask", landIceMask)

          call MPAS_pool_get_array(meshPool, "latCell", latCell)
          call MPAS_pool_get_array(meshPool, "lonCell", lonCell)

          call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

          do iCell = 1, nCells

             if (lonCell(iCell) > lonLimitWest * ciceDegreesToRadians .and. &
                 lonCell(iCell) < lonLimitEast * ciceDegreesToRadians .and. &
                 latCell(iCell) < latLimit     * ciceDegreesToRadians) then

                landIceMask(iCell) = 1

             else

                landIceMask(iCell) = 0

             endif

          enddo ! iCell

          block => block % next
       enddo

    endif

  end subroutine init_test_ice_shelf_mask

!-----------------------------------------------------------------------

end module cice_initialize
