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

  implicit none

  private
  save

  public :: &
       cice_init

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

    use cice_mesh, only: cice_init_boundary
    use cice_advection, only: cice_init_advection
    use cice_velocity_solver, only: cice_init_velocity_solver
    use cice_debug, only: cice_test_case_output, cice_debug_output_domain
    use cice_column, only: cice_init_column_physics_package

    type(domain_type), intent(inout) :: &
         domain !< Input/Output: 

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input: 

    real(kind=RKIND), intent(in) :: &
         dt !< Input: 

    type(block_type), pointer :: &
         block

    logical, pointer :: &
         config_use_velocity_solver, &
         config_use_advection, &
         config_use_forcing, &
         config_do_restart, &
         config_use_column_physics_package

    call MPAS_pool_get_config(domain % configs, "config_use_velocity_solver", config_use_velocity_solver)
    call MPAS_pool_get_config(domain % configs, "config_use_advection", config_use_advection)
    call MPAS_pool_get_config(domain % configs, "config_use_forcing", config_use_forcing)
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)
    call MPAS_pool_get_config(domain % configs, "config_use_column_physics_package", config_use_column_physics_package)

    block => domain % blocklist
    do while (associated(block))
       call init_junk_values(block)
       block => block % next
    end do
    
    ! calculate boundary variables
    call cice_init_boundary(domain)

    ! init dynamics
    if (config_use_velocity_solver) then
       block => domain % blocklist
       do while (associated(block))
          call cice_init_velocity_solver(block, dt)
          block => block % next
       end do
    endif
    
    ! init advection
    if (config_use_advection) then
       block => domain % blocklist
       do while (associated(block))
          call cice_init_advection(block)
          block => block % next
       end do
    endif

    ! column physics initialization
    if (config_use_column_physics_package) then
       call cice_init_column_physics_package(&
            domain, &
            clock)
    endif

    ! init ice state
    if (.not. config_do_restart) &
         call init_ice_state(domain)

    call cice_debug_output_domain(domain)

    block => domain % blocklist
    do while (associated(block))
       call cice_test_case_output(block, 0)
       block => block % next
    end do

    ! initial halo exchange
    call initial_halo_exchanges(domain)

  end subroutine cice_init!}}}

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

  subroutine init_junk_values(block)!{{{

    use cice_debug, only: &
         cice_init_ijpop_from_ivertex

    type(block_type), intent(inout) :: &
         block !< Input/Output: 

    type(MPAS_pool_type), pointer :: &
         mesh

    integer :: &
         njunk

    real(kind=RKIND), parameter :: &
         value_junk = -1.0e34

    integer, pointer :: &
         nCells
    
    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

    !call cice_init_ijpop_from_ivertex(mesh)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    njunk = nCells + 1
    
    areaCell(njunk) = value_junk

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
       domain)!{{{

    use cice_testing, only: &
         cice_init_square_test_case_hex, &
         cice_init_square_point_test_case_hex

    type(domain_type), intent(inout) :: &
         domain !< Input/Output: 

    type(block_type), pointer :: &
         block

    character(len=strKIND), pointer :: &
         config_initial_condition_type, &
         config_initial_velocity_type

    integer, dimension(:), pointer :: &
         interiorVertex

    type(MPAS_pool_type), pointer :: &
         tracers_fields, &
         velocity_solver_fields

    type(field3DReal), pointer :: &
         iceAreaCategoryField, &
         iceVolumeCategoryField, &
         snowVolumeCategoryField, &
         surfaceTemperatureField

    type(field1DReal), pointer :: &
         uVelocityField, &
         vVelocityField

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

       else

          write(*,*) "Ice initial condition type unknown!"
          write(*,*) "config_initial_condition_type: ", trim(config_initial_condition_type)
          stop

       endif

       block => block % next
    end do

    ! halo exchanges
    call MPAS_pool_get_subpool(domain % blocklist % structs, "tracers", tracers_fields)
    
    call MPAS_pool_get_field(tracers_fields, "iceAreaCategory", iceAreaCategoryField)
    call MPAS_pool_get_field(tracers_fields, "iceVolumeCategory", iceVolumeCategoryField)
    call MPAS_pool_get_field(tracers_fields, "snowVolumeCategory", snowVolumeCategoryField)
    call MPAS_pool_get_field(tracers_fields, "surfaceTemperature", surfaceTemperatureField)

    call MPAS_dmpar_exch_halo_field(iceAreaCategoryField)
    call MPAS_dmpar_exch_halo_field(iceVolumeCategoryField)
    call MPAS_dmpar_exch_halo_field(snowVolumeCategoryField)
    call MPAS_dmpar_exch_halo_field(surfaceTemperatureField)

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

       else

          write(*,*) "Ice initial velocity type unknown!"
          write(*,*) "config_initial_velocity_type: ", trim(config_initial_velocity_type)
          stop

       endif

       block => block % next
    end do

    ! halo exchanges
    call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocity_solver_fields)

    call MPAS_pool_get_field(velocity_solver_fields, "uVelocity", uVelocityField)
    call MPAS_pool_get_field(velocity_solver_fields, "vVelocity", vVelocityField)

    call MPAS_dmpar_exch_halo_field(uVelocityField)
    call MPAS_dmpar_exch_halo_field(vVelocityField)

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
         icestate

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
    call MPAS_pool_get_subpool(block % structs, "icestate", icestate)

    call MPAS_pool_get_dimension(mesh, "nCells", nCellsSolve)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    call MPAS_pool_get_array(icestate, "iceAreaCell", iceAreaCell)
    call MPAS_pool_get_array(icestate, "surfaceTemperatureCell", surfaceTemperatureCell)

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
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature)
    call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy)
    call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy)
    call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity, 1)

    call MPAS_pool_get_array(atmos_coupling, "airTemperature", airTemperature)

    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTemperature", seaSurfaceTemperature)
    call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)

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
           latCell(iCell) < config_initial_latitude_south * ciceDegreesToRadians)) then

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

    ! clean up
    deallocate(initialCategoryIceArea)
    deallocate(initialCategoryIceThickness)

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

    real(kind=RKIND), parameter :: &
         thicknessWithLargestArea = 3.0_RKIND ! initial ice thickness with greatest area

    call MPAS_pool_get_subpool(block % structs, "initial", initial)
    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

    call MPAS_pool_get_array(initial, "categoryThicknessLimits", categoryThicknessLimits)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

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

    integer :: &
         iCell
    
    real(kind=RKIND), parameter :: &
         circle_radius = 1.0e6_RKIND
   
    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)

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

    do iCell = 1, nCells
       
       if (sqrt(xCell(iCell)**2 + zCell(iCell)**2) < circle_radius .and. yCell(iCell) > 0.0_RKIND) then
          
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

    do iVertex = 1, nVertices

       if (interiorVertex(iVertex) == 1) then
          
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
         ocean_forcing

    integer, pointer :: &
         nVertices

    integer, dimension(:), pointer :: &
         interiorVertex

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
    call MPAS_pool_get_subpool(block % structs, "ocean_forcing", ocean_forcing)

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)

    call MPAS_pool_get_array(ocean_forcing, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(ocean_forcing, "vOceanVelocity", vOceanVelocity)

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
    
       if (interiorVertex(iVertex) == 1) then

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

    type(MPAS_pool_type), pointer :: &
         velocity_weak_fields

    type(field3DReal), pointer :: &
         normalVectorTriangleField

    call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_weak", velocity_weak_fields)

    call MPAS_pool_get_field(velocity_weak_fields, "normalVectorTriangle", normalVectorTriangleField)

    call MPAS_dmpar_exch_halo_field(normalVectorTriangleField)

  end subroutine initial_halo_exchanges

!-----------------------------------------------------------------------

end module cice_initialize
