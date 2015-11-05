!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_velocity_solver
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_velocity_solver

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_dmpar
  use mpas_timer

  implicit none

  private
  save

  public :: &
       cice_init_velocity_solver, &
       cice_run_velocity_solver

  ! velocity solver constants
  real(kind=RKIND), parameter, private :: &
       sinOceanTurningAngle = 0.0_RKIND, & ! northern hemisphere
       cosOceanTurningAngle = 1.0_RKIND, & ! northern hemisphere
       ciceAreaMinimum = 0.001_RKIND, &
       ciceMassMinimum = 0.01_RKIND

contains

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_velocity_solver
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_init_velocity_solver(&
       domain)!{{{

    use cice_velocity_solver_constitutive_relation, only: &
         cice_init_evp

    use cice_velocity_solver_weak, only: &
         cice_init_velocity_solver_weak

    use cice_velocity_solver_variational, only: &
         cice_init_velocity_solver_variational

    type (domain_type), intent(inout) :: &
         domain !< Input/Output: 

    type(block_type), pointer :: &
         block

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme, &
         config_variational_basis

    logical, pointer :: &
         config_use_velocity_solver, &
         config_rotate_cartesian_grid, &
         config_include_metric_terms

    type (MPAS_pool_type), pointer :: &
         mesh, &
         boundary, &
         velocitySolver, &
         velocity_weak, &
         velocity_variational, &
         velocity_pwl

    real(kind=RKIND), pointer :: &
         dynamicsTimeStep, &
         elasticTimeStep, &
         config_dt

    integer, pointer :: &
         config_dynamics_subcycle_number, &
         config_elastic_subcycle_number

    ! set timesteps even with velocity turned off
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_dynamics_subcycle_number", config_dynamics_subcycle_number)
       call MPAS_pool_get_config(block % configs, "config_elastic_subcycle_number", config_elastic_subcycle_number)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolver)
       call MPAS_pool_get_array(velocitySolver, "dynamicsTimeStep", dynamicsTimeStep)
       call MPAS_pool_get_array(velocitySolver, "elasticTimeStep", elasticTimeStep)

       dynamicsTimeStep = config_dt / real(config_dynamics_subcycle_number,RKIND)

       elasticTimeStep = dynamicsTimeStep / real(config_elastic_subcycle_number,RKIND)

       block => block % next
    enddo

    ! check if we initialize velocity solver
    call MPAS_pool_get_config(domain % configs, "config_use_velocity_solver", config_use_velocity_solver)

    if (config_use_velocity_solver) then
 
       ! initialize the evp solver
       call cice_init_evp(domain)

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
          call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
          call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocity_weak)
          call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocity_variational)

          call MPAS_pool_get_config(block % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)
          call MPAS_pool_get_config(block % configs, "config_variational_basis", config_variational_basis)
          call MPAS_pool_get_config(block % configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)
          call MPAS_pool_get_config(block % configs, "config_include_metric_terms", config_include_metric_terms)

          if (trim(config_stress_divergence_scheme) == "weak") then

             call cice_init_velocity_solver_weak(&
                  mesh, &
                  boundary, &
                  velocity_weak, &
                  config_rotate_cartesian_grid)

          else if (trim(config_stress_divergence_scheme) == "variational") then

             call cice_init_velocity_solver_variational(&
                  mesh, &
                  velocity_variational, &
                  boundary, &
                  config_rotate_cartesian_grid, &
                  config_include_metric_terms, &
                  config_variational_basis)

          else if (trim(config_stress_divergence_scheme) == "none") then

             continue

          else 

             write(stderrUnit,*) "Error: config_stress_divergence_scheme unknown: ", trim(config_stress_divergence_scheme)
             stop

          endif

          block => block % next
       enddo

    endif

  end subroutine cice_init_velocity_solver!}}}

!-----------------------------------------------------------------------
! Time stepping
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_run_velocity_solver
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_run_velocity_solver(domain, clock)!{{{

    use cice_debug

    type (domain_type), intent(inout) :: &
         domain !< Input/Output: 

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input: 

    logical, pointer :: &
         config_use_velocity_solver

    ! determine if velocity solver switched on
    call MPAS_pool_get_config(domain % configs, "config_use_velocity_solver", config_use_velocity_solver)

    if (config_use_velocity_solver) then

       ! pre subcycle
       call mpas_timer_start("Velocity solver pre-cycle")
       call velocity_solver_pre_subcycle(domain)
       call mpas_timer_stop("Velocity solver pre-cycle")

       ! subcycle the dynamics
       call mpas_timer_start("Velocity solver sub-cycle")
       call subcycle_velocity_solver(domain, clock)
       call mpas_timer_stop("Velocity solver sub-cycle")

       ! post subcycle
       call mpas_timer_start("Velocity solver post-cycle")
       call velocity_solver_post_subcycle(domain)
       call mpas_timer_stop("Velocity solver post-cycle")

    endif ! config_use_velocity_solver

  end subroutine cice_run_velocity_solver!}}}

!-----------------------------------------------------------------------
! Pre sub-cycle
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  velocity_solver_pre_subcycle
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date January 13th 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine velocity_solver_pre_subcycle(domain)

    type(domain_type), intent(inout) :: &
         domain

    ! aggregate categories for area and volume into total mass
    call aggregate_mass_and_area(domain)

    ! calculate computational masks
    call calculation_masks(domain)

    ! set new ice velocities to ocean velocity
    call new_ice_velocities(domain)

    ! calculate the ice strength
    call ice_strength(domain)

    ! calculate air stress
    call air_stress(domain)

    ! calculate the coriolis force coefficient
    call coriolis_force_coefficient(domain)

    ! calculate the ocean stress
    call ocean_stress(domain)

    ! calculate the surface tilt force
    call surface_tilt(domain)

    ! initialize subcycle variables
    call init_subcycle_variables(domain)

  end subroutine velocity_solver_pre_subcycle

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  aggregate_mass_and_area
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine aggregate_mass_and_area(domain)!{{{

    use cice_constants, only: &
         ciceDensityIce, &
         ciceDensitySnow

    type(domain_type), intent(in) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         tracersPool, &
         tracersAggregatePool, &
         icestatePool

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    real(kind=RKIND), dimension(:), pointer :: &
         totalMassCell, &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell

    type(field3DReal), pointer :: &
         iceAreaCategoryField, &
         iceVolumeCategoryField, &
         snowVolumeCategoryField

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    ! halo exchange of mass tracers
    call mpas_timer_start("Halo")

    call MPAS_pool_get_subpool(domain % blocklist % structs, "tracers", tracersPool)

    call MPAS_pool_get_field(tracersPool, "iceAreaCategory", iceAreaCategoryField, 1)
    call MPAS_pool_get_field(tracersPool, "iceVolumeCategory", iceVolumeCategoryField, 1)
    call MPAS_pool_get_field(tracersPool, "snowVolumeCategory", snowVolumeCategoryField, 1)

    call MPAS_dmpar_exch_halo_field(iceAreaCategoryField)
    call MPAS_dmpar_exch_halo_field(iceVolumeCategoryField)
    call MPAS_dmpar_exch_halo_field(snowVolumeCategoryField)

    call mpas_timer_stop("Halo")

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

       call MPAS_pool_get_subpool(block % structs, "tracers", tracersPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

       call MPAS_pool_get_array(tracersPool, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracersPool, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracersPool, "snowVolumeCategory", snowVolumeCategory, 1)

       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(tracersAggregatePool, "iceVolumeCell", iceVolumeCell)
       call MPAS_pool_get_array(tracersAggregatePool, "snowVolumeCell", snowVolumeCell)

       call MPAS_pool_get_array(icestatePool, "totalMassCell", totalMassCell)

       do iCell = 1, nCells

          iceAreaCell(iCell)    = sum(iceAreaCategory(1,:,iCell))
          iceVolumeCell(iCell)  = sum(iceVolumeCategory(1,:,iCell))
          snowVolumeCell(iCell) = sum(snowVolumeCategory(1,:,iCell))

          totalMassCell(iCell)  = iceVolumeCell(iCell)  * ciceDensityIce + &
                                  snowVolumeCell(iCell) * ciceDensitySnow

       enddo ! iCell

       block => block % next
    enddo

  end subroutine aggregate_mass_and_area!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calculation_masks
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 19th September 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine calculation_masks(domain)

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         icestatePool, &
         tracersAggregatePool, &
         meshPool, &
         velocitySolverPool

    type(field1DReal), pointer :: &
         iceAreaCellInitialField

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaVertex, &
         iceAreaCellInitial, &
         iceAreaCell, &
         totalMassVertex, &
         totalMassCell

    type(field1DInteger), pointer :: &
         solveVelocityField

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_vertical_thermodynamics

    ! set initial ice area if not have column physics
    call MPAS_pool_get_config(domain % configs, "config_use_column_package", config_use_column_package)
    call MPAS_pool_get_config(domain % configs, "config_use_column_vertical_thermodynamics", config_use_column_vertical_thermodynamics)

    if (.not. config_use_column_package .or. &
         (config_use_column_package .and. .not. config_use_column_vertical_thermodynamics)) then

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
          call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
          call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
          call MPAS_pool_get_array(icestatePool, "iceAreaCellInitial", iceAreaCellInitial)

          iceAreaCellInitial = iceAreaCell

          block => block % next
       end do

    endif

    ! halo exchange of initial ice area of cell
    call mpas_timer_start("Halo")
    call MPAS_pool_get_subpool(domain % blocklist % structs, "icestate", icestatePool)
    call MPAS_pool_get_field(icestatePool, "iceAreaCellInitial", iceAreaCellInitialField)
    call MPAS_dmpar_exch_halo_field(iceAreaCellInitialField)
    call mpas_timer_stop("Halo")

    ! interpolate area and mass from cells to vertices
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

       call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)
       call MPAS_pool_get_array(icestatePool, "totalMassCell", totalMassCell)
       call MPAS_pool_get_array(icestatePool, "iceAreaCellInitial", iceAreaCellInitial)
       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            iceAreaVertex, &
            iceAreaCellInitial)

       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            totalMassVertex, &
            totalMassCell)

       block => block % next
    end do

    ! calculate computational masks
    call stress_calculation_mask(domain)

    call velocity_calculation_mask(domain)

    ! halo exchange velocity mask
    call mpas_timer_start("Halo")
    call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolverPool)
    call MPAS_pool_get_field(velocitySolverPool, "solveVelocity", solveVelocityField)
    call MPAS_dmpar_exch_halo_field(solveVelocityField)
    call mpas_timer_stop("Halo")

  end subroutine calculation_masks

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  stress_calculation_mask
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine stress_calculation_mask(domain)!{{{

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocitySolverPool, &
         icestatePool

    integer, dimension(:), pointer :: &
         solveStress

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCellInitial, &
         totalMassCell

    integer :: &
         iCell, &
         iCellOnCell, &
         iCellNeighbour

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell
    
    integer, dimension(:,:), pointer :: &
         cellsOnCell

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

       call MPAS_pool_get_array(meshPool, "nEdgesOnCell", nEdgesOnCell)
       call MPAS_pool_get_array(meshPool, "cellsOnCell", cellsOnCell)

       call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)

       call MPAS_pool_get_array(icestatePool, "iceAreaCellInitial", iceAreaCellInitial)
       call MPAS_pool_get_array(icestatePool, "totalMassCell", totalMassCell)

       do iCell = 1, nCells

          solveStress(iCell) = 0

          if (iceAreaCellInitial(iCell) > ciceAreaMinimum .and. &
              totalMassCell(iCell) > ciceMassMinimum) then

             ! this cell has sufficient ice
             solveStress(iCell) = 1

          else

             ! test neighbouring cells to see if have sufficient ice
             do iCellOnCell = 1, nEdgesOnCell(iCell)

                iCellNeighbour = cellsOnCell(iCellOnCell,iCell)

                if (iceAreaCellInitial(iCellNeighbour) > ciceAreaMinimum .and. &
                    totalMassCell(iCellNeighbour) > ciceMassMinimum) then

                   solveStress(iCell) = 1
                   exit

                endif

             enddo ! iCellOnCell

          endif

       enddo ! iCell

       do iCell = nCells+1, nCells

          solveStress(iCell) = 0

       enddo ! iCell

       block => block % next
    enddo

  end subroutine stress_calculation_mask!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  velocity_calculation_mask
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine velocity_calculation_mask(domain)!{{{

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         icestatePool, &
         boundaryPool

    integer, dimension(:), pointer :: &
         solveVelocity

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaVertex, &
         totalMassVertex

    integer, dimension(:), pointer :: &
         interiorVertex

    integer :: &
         iVertex

    integer, pointer :: &
         nVerticesSolve, &
         nVertices

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nVertices", nVertices)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "boundary", boundaryPool)

       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)

       call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)
       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

       call MPAS_pool_get_array(boundaryPool, "interiorVertex", interiorVertex)

       do iVertex = 1, nVerticesSolve

          solveVelocity(iVertex) = 0

          if (interiorVertex(iVertex) == 1 .and. &
               iceAreaVertex(iVertex) > ciceAreaMinimum .and. &
               totalMassVertex(iVertex) > ciceMassMinimum) then

             ! this vertex has sufficient ice
             solveVelocity(iVertex) = 1

          endif

       enddo ! iVertex

       do iVertex = nVerticesSolve+1, nVertices

          solveVelocity(iVertex) = 0

       enddo ! iVertex

       block => block % next
    enddo

  end subroutine velocity_calculation_mask!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  new_ice_velocities
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 29th June 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine new_ice_velocities(domain)!{{{

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         meshPool, &
         oceanCouplingPool

    real(kind=RKIND), dimension(:), pointer :: &
         uOceanVelocity, &
         vOceanVelocity, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         uVelocity, &
         vVelocity, &
         uVelocityInitial, &
         vVelocityInitial, &
         stressDivergenceU, &
         stressDivergenceV, &
         oceanStressU, &
         oceanStressV

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveVelocityPrevious

    type(field1DReal), pointer :: &
         uVelocityField, &
         vVelocityField

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)

       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocityInitial", uVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocityInitial", vVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocityPrevious", solveVelocityPrevious)
       call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
       call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)

       call MPAS_pool_get_array(oceanCouplingPool, "uOceanVelocity", uOceanVelocity)
       call MPAS_pool_get_array(oceanCouplingPool, "vOceanVelocity", vOceanVelocity)

       ! interpolate cell ocean velocity to vertex
       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            uOceanVelocityVertex, &
            uOceanVelocity)

       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            vOceanVelocityVertex, &
            vOceanVelocity)

       ! set new ice to ocean velocity
       do iVertex = 1, nVerticesSolve

          if (solveVelocity(iVertex) == 1) then

             if (solveVelocityPrevious(iVertex) == 0) then

                uVelocity(iVertex) = uOceanVelocityVertex(iVertex)
                vVelocity(iVertex) = vOceanVelocityVertex(iVertex)

             endif

          else

             uVelocity(iVertex) = 0.0_RKIND
             vVelocity(iVertex) = 0.0_RKIND
             stressDivergenceU(iVertex) = 0.0_RKIND
             stressDivergenceV(iVertex) = 0.0_RKIND
             oceanStressU(iVertex) = 0.0_RKIND
             oceanStressV(iVertex) = 0.0_RKIND

          endif

       enddo ! iVertex

       solveVelocityPrevious = solveVelocity

       uVelocityInitial = uVelocity
       vVelocityInitial = vVelocity

       block => block % next
    enddo

    ! halo exchange velocities
    call mpas_timer_start("Halo")

    call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolverPool)
    call MPAS_pool_get_field(velocitySolverPool, "uVelocity", uVelocityField)
    call MPAS_pool_get_field(velocitySolverPool, "vVelocity", vVelocityField)
    call MPAS_dmpar_exch_halo_field(uVelocityField)
    call MPAS_dmpar_exch_halo_field(vVelocityField)

    call mpas_timer_stop("Halo")

  end subroutine new_ice_velocities!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ice_strength
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine ice_strength(domain)

    use ice_colpkg, only: &
         colpkg_ice_strength

    use cice_constants, only: &
         ciceIceStrengthConstantHiblerP, &
         ciceIceStrengthConstantHiblerC

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         tracersAggregatePool, &
         icestatePool, &
         tracersPool

    real(kind=RKIND), dimension(:), pointer:: &
         icePressure, &
         iceAreaCell, &
         iceVolumeCell, &
         openWaterArea

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory

    integer, dimension(:), pointer :: &
         solveStress

    type(field1DReal), pointer :: &
         icePressureField

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_vertical_thermodynamics

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    integer :: &
         iCell

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_column_package", config_use_column_package)
       call MPAS_pool_get_config(block % configs, "config_use_column_vertical_thermodynamics", config_use_column_vertical_thermodynamics)

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nCategories", nCategories)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracersPool)

       call MPAS_pool_get_array(velocitySolverPool, "icePressure", icePressure)
       call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)

       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(tracersAggregatePool, "iceVolumeCell", iceVolumeCell)

       call MPAS_pool_get_array(icestatePool, "openWaterArea", openWaterArea)

       call MPAS_pool_get_array(tracersPool, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracersPool, "iceVolumeCategory", iceVolumeCategory, 1)

       if (.not. config_use_column_package .or. &
            (config_use_column_package .and. .not. config_use_column_vertical_thermodynamics)) then

          do iCell = 1, nCellsSolve

             if (solveStress(iCell) == 1) then

                icePressure(iCell) = ciceIceStrengthConstantHiblerP * iceVolumeCell(iCell) * &
                     exp(-ciceIceStrengthConstantHiblerC*(1.0_RKIND-iceAreaCell(iCell)))

             else

                icePressure(iCell) = 0.0_RKIND

             endif ! solveStress

          enddo ! iCell

       else

          do iCell = 1, nCellsSolve

             icePressure(iCell) = 0.0_RKIND

             if (solveStress(iCell) == 1) then

                ! this routine doesnt reset icePressure
                call colpkg_ice_strength(&
                     nCategories, &
                     iceAreaCell(iCell), &
                     iceVolumeCell(iCell), &
                     openWaterArea(iCell), &
                     iceAreaCategory(1,:,iCell), &
                     iceVolumeCategory(1,:,iCell), &
                     icePressure(iCell))

             endif ! solveStress

          enddo ! iCell

       endif

       block => block % next
    enddo

    ! halo exchange ice strength
    call mpas_timer_start("Halo")

    call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolverPool)
    call MPAS_pool_get_field(velocitySolverPool, "icePressure", icePressureField)
    call MPAS_dmpar_exch_halo_field(icePressureField)

    call mpas_timer_stop("Halo")

  end subroutine ice_strength

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  air_stress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 19th September 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine air_stress(domain)

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_vertical_thermodynamics, &
         config_use_air_stress

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocitySolverPool

    real(kind=RKIND), dimension(:), pointer :: &
         airStressVertexU, &
         airStressVertexV, &
         airStressCellU, &
         airStressCellV

    type(field1DReal), pointer :: &
         airStressCellUField, &
         airStressCellVField!, &

    ! calculate the air stress
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_package", config_use_column_package)
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_vertical_thermodynamics", config_use_column_vertical_thermodynamics)

    if (.not. config_use_column_package .or. &
         (config_use_column_package .and. .not. config_use_column_vertical_thermodynamics)) then
       call constant_air_stress(domain)
    endif

    ! check for no air stress
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_air_stress", config_use_air_stress)
    if (.not. config_use_air_stress) then
       block => domain % blocklist
       do while (associated(block))
          call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
          call MPAS_pool_get_array(velocitySolverPool, "airStressCellU", airStressCellU)
          call MPAS_pool_get_array(velocitySolverPool, "airStressCellV", airStressCellV)
          airStressCellU = 0.0_RKIND
          airStressCellV = 0.0_RKIND
          block => block % next
       end do
    endif ! .not. config_use_air_stress

    ! halo exchange air stress
    call mpas_timer_start("Halo")

    call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolverPool)
    call MPAS_pool_get_field(velocitySolverPool, "airStressCellU", airStressCellUField)
    call MPAS_pool_get_field(velocitySolverPool, "airStressCellV", airStressCellVField)
    call MPAS_dmpar_exch_halo_field(airStressCellUField)
    call MPAS_dmpar_exch_halo_field(airStressCellVField)

    call mpas_timer_stop("Halo")

    ! interpolate air stress
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressCellV", airStressCellV)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexU", airStressVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexV", airStressVertexV)

       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            airStressVertexU, &
            airStressCellU)

       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            airStressVertexV, &
            airStressCellV)

       block => block % next
    end do

  end subroutine air_stress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  constant_air_stress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine constant_air_stress(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         atmosCouplingPool, &
         tracersAggregatePool

    real(kind=RKIND), dimension(:), pointer :: &
         airStressCellU, &
         airStressCellV, &
         uAirVelocity, &
         vAirVelocity, &
         airDensity, &
         iceAreaCell

    real(kind=RKIND) :: &
         windSpeed

    integer, pointer :: &
         nCellsSolve

    integer :: &
         iCell

    real(kind=RKIND), parameter :: &
         airStressCoeff = 0.0012_RKIND 

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCouplingPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)

       call MPAS_pool_get_array(velocitySolverPool, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressCellV", airStressCellV)

       call MPAS_pool_get_array(atmosCouplingPool, "uAirVelocity", uAirVelocity)
       call MPAS_pool_get_array(atmosCouplingPool, "vAirVelocity", vAirVelocity)
       call MPAS_pool_get_array(atmosCouplingPool, "airDensity", airDensity)

       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)

       do iCell = 1, nCellsSolve

          windSpeed = sqrt(uAirVelocity(iCell)**2 + vAirVelocity(iCell)**2)

          airStressCellU(iCell) = airDensity(iCell) * windSpeed * airStressCoeff * uAirVelocity(iCell) * iceAreaCell(iCell)
          airStressCellV(iCell) = airDensity(iCell) * windSpeed * airStressCoeff * vAirVelocity(iCell) * iceAreaCell(iCell)

       enddo ! iCell

       block => block % next
    end do

  end subroutine constant_air_stress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  coriolis_force_coefficient
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine coriolis_force_coefficient(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         icestatePool, &
         velocitySolverPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         totalMassVertexfVertex, &
         totalMassVertex, &
         fVertex

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)
       call MPAS_pool_get_array(velocitySolverPool, "totalMassVertexfVertex", totalMassVertexfVertex)
       call MPAS_pool_get_array(meshPool, "fVertex", fVertex)

       do iVertex = 1, nVerticesSolve

          totalMassVertexfVertex(iVertex) = totalMassVertex(iVertex) * fVertex(iVertex)

       enddo ! iVertex

       block => block % next
    end do

  end subroutine coriolis_force_coefficient

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocean_stress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine ocean_stress(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         oceanStressU, &
         oceanStressV, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         fVertex

    integer, dimension(:), pointer :: &
         solveVelocity

    logical, pointer :: &
         configUseOceanStress

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_ocean_stress", configUseOceanStress)

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)

       if (configUseOceanStress) then

          call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
          call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)

          call MPAS_pool_get_array(meshPool, "fVertex", fVertex)

          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) == 1) then

                oceanStressU(iVertex) = uOceanVelocityVertex(iVertex) * cosOceanTurningAngle - &
                                        vOceanVelocityVertex(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,fVertex(iVertex))
                oceanStressV(iVertex) = uOceanVelocityVertex(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,fVertex(iVertex)) + &
                                        vOceanVelocityVertex(iVertex) * cosOceanTurningAngle

             else

                oceanStressU(iVertex) = 0.0_RKIND
                oceanStressV(iVertex) = 0.0_RKIND

             endif ! solvePoints

          enddo ! iVertex

       else

          ! no ocean stress
          oceanStressU = 0.0_RKIND
          oceanStressV = 0.0_RKIND

       endif

       block => block % next
    end do

  end subroutine ocean_stress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  surface_tilt
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine surface_tilt(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         icestatePool, &
         velocitySolverPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer:: &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         totalMassVertex, &
         fVertex

    integer, dimension(:), pointer :: &
         solveVelocity

    logical, pointer :: &
         configUseSurfaceTilt

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_surface_tilt", configUseSurfaceTilt)

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)

       if (configUseSurfaceTilt) then

          call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
          call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)

          call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

          call MPAS_pool_get_array(meshPool, "fVertex", fVertex)

          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) == 1) then

                surfaceTiltForceU(iVertex) = -fVertex(iVertex) * totalMassVertex(iVertex) * vOceanVelocityVertex(iVertex)
                surfaceTiltForceV(iVertex) =  fVertex(iVertex) * totalMassVertex(iVertex) * uOceanVelocityVertex(iVertex)

             else

                surfaceTiltForceU(iVertex) = 0.0_RKIND
                surfaceTiltForceV(iVertex) = 0.0_RKIND

             endif ! solveVelocity

          enddo ! iVertex

       else

          ! no surface tilt
          surfaceTiltForceU = 0.0_RKIND
          surfaceTiltForceV = 0.0_RKIND

       endif

       block => block % next
    end do

  end subroutine surface_tilt

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_subcycle_variables
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 18th October 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_subcycle_variables(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         velocityVariationalPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11, &
         strain22, &
         strain12, &
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND), dimension(:), pointer :: &
         stressDivergenceU, &
         stressDivergenceV, &
         uVelocity, &
         vVelocity, &
         oceanStressCoeff

    integer, dimension(:), pointer :: &
         solveStress, &
         solveVelocity

    integer, pointer :: &
         nCells, &
         nVerticesSolve

    integer :: &
         iCell, &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocityVariationalPool)

       ! strains
       call MPAS_pool_get_array(velocityVariationalPool, "strain11", strain11)
       call MPAS_pool_get_array(velocityVariationalPool, "strain22", strain22)
       call MPAS_pool_get_array(velocityVariationalPool, "strain12", strain12)

       strain11 = 0.0_RKIND
       strain22 = 0.0_RKIND
       strain12 = 0.0_RKIND

       ! stresses
       call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)

       call MPAS_pool_get_array(velocityVariationalPool, "stress11", stress11)
       call MPAS_pool_get_array(velocityVariationalPool, "stress22", stress22)
       call MPAS_pool_get_array(velocityVariationalPool, "stress12", stress12)

       call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

       do iCell = 1, nCells

          if (solveStress(iCell) /= 1) then

             stress11(:,iCell) = 0.0_RKIND
             stress22(:,iCell) = 0.0_RKIND
             stress12(:,iCell) = 0.0_RKIND

          endif

       enddo ! iCell

       ! divergence of stress
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)

       stressDivergenceU = 0.0_RKIND
       stressDivergenceV = 0.0_RKIND

       ! sea ice velocity
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)

       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       do iVertex = 1, nVerticesSolve

          if (solveVelocity(iVertex) /= 1) then

             uVelocity(iVertex) = 0.0_RKIND
             vVelocity(iVertex) = 0.0_RKIND

             oceanStressCoeff(iVertex) = 0.0_RKIND

          endif

       enddo ! iVertex

       block => block % next
    enddo

  end subroutine init_subcycle_variables

!-----------------------------------------------------------------------
! Sub-cycle
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  subcycle_velocity_solver
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine subcycle_velocity_solver(&
       domain, &
       clock)!{{{

    type(domain_type), intent(inout) :: &
         domain !< Input/Output: 

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input: 
 
    integer, pointer :: &
         config_elastic_subcycle_number

    integer :: &
         iElasticSubcycle

    call MPAS_pool_get_config(domain % configs, "config_elastic_subcycle_number", config_elastic_subcycle_number)

    do iElasticSubcycle = 1, config_elastic_subcycle_number

       call single_subcycle_velocity_solver(&
            domain, &
            clock, &
            iElasticSubcycle)

    enddo

  end subroutine subcycle_velocity_solver!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  single_subcycle_velocity_solver
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine single_subcycle_velocity_solver(&
       domain, &
       clock, &
       iElasticSubcycle)!{{{

    use cice_velocity_solver_weak, only: &
         cice_internal_stress_weak

    use cice_velocity_solver_variational, only: &
         cice_internal_stress_variational

    use cice_debug, only: &
         cice_time_output_variable_real

    type(domain_type), intent(inout) :: &
         domain !< Input/Output: 

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input: 

    integer, intent(in) :: &
         iElasticSubcycle !< Input: !! testing

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    logical, pointer :: &
         config_revised_evp

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool

    type(field1DReal), pointer :: &
         uVelocity, &
         vvelocity

    call MPAS_pool_get_config(domain % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)
    call MPAS_pool_get_config(domain % configs, "config_revised_evp", config_revised_evp)

    ! calculate internal stresses
    if (trim(config_stress_divergence_scheme) == "weak") then

       call mpas_timer_start("Velocity solver internal stress")
       call cice_internal_stress_weak(domain)
       call mpas_timer_stop("Velocity solver internal stress")

    else if (trim(config_stress_divergence_scheme) == "variational") then

       call mpas_timer_start("Velocity solver internal stress")
       call cice_internal_stress_variational(domain)
       call mpas_timer_stop("Velocity solver internal stress")

    endif

    ! ocean stress coefficient
    call ocean_stress_coefficient(domain)

    ! solve for velocity
    if (.not. config_revised_evp) then

       call mpas_timer_start("Velocity solver compute")
       call solve_velocity(domain)
       call mpas_timer_stop("Velocity solver compute")

    else

       call mpas_timer_start("Velocity solver compute")
       call solve_velocity_revised(domain)
       call mpas_timer_stop("Velocity solver compute")

    endif

    ! halo exchange
    call mpas_timer_start("Halo")
    call mpas_timer_start("Velocity solver halo")

    call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolverPool)

    call MPAS_pool_get_field(velocitySolverPool, "uVelocity", uVelocity)
    call MPAS_pool_get_field(velocitySolverPool, "vVelocity", vVelocity)

    call MPAS_dmpar_exch_halo_field(uVelocity)
    call MPAS_dmpar_exch_halo_field(vVelocity)

    call mpas_timer_stop("Velocity solver halo")
    call mpas_timer_stop("Halo")

  end subroutine single_subcycle_velocity_solver!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocean_stress_coefficient
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine ocean_stress_coefficient(domain)

    use cice_constants, only: &
         ciceDensitySeaWater, &
         ciceIceOceanDragCoefficient

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         icestatePool

    real(kind=RKIND), dimension(:), pointer :: &
         oceanStressCoeff, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         uVelocity, &
         vVelocity, &
         iceAreaVertex

    integer, dimension(:), pointer :: &
         solveVelocity

    logical, pointer :: &
         configUseOceanStress

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_ocean_stress", configUseOceanStress)

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)

       if (configUseOceanStress) then

          call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
          call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)
          call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)

          call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)

          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) == 1) then

                oceanStressCoeff(iVertex) = ciceIceOceanDragCoefficient * ciceDensitySeaWater * iceAreaVertex(iVertex) * &
                     sqrt((uOceanVelocityVertex(iVertex) - uVelocity(iVertex))**2 + &
                          (vOceanVelocityVertex(iVertex) - vVelocity(iVertex))**2)

             endif

          enddo ! iVertex

       else

          ! no ocean stress
          oceanStressCoeff = 0.0_RKIND

       endif

       block => block % next
    end do

  end subroutine ocean_stress_coefficient

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  solve_velocity
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine solve_velocity(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         icestatePool

    integer, dimension(:), pointer :: &
         solveVelocity

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         uVelocityInitial, &
         vVelocityInitial, &
         totalMassVertex, &
         totalMassVertexfVertex, &
         stressDivergenceU, &
         stressDivergenceV, &
         airStressVertexU, &
         airStressVertexV, &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         oceanStressU, &
         oceanStressV, &
         oceanStressCoeff

    real(kind=RKIND), pointer :: &
         elasticTimeStep

    real(kind=RKIND), dimension(2) :: &
         rightHandSide

    real(kind=RKIND), dimension(2,2) :: &
         leftMatrix

    real(kind=RKIND) :: &
         solutionDenominator

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

       call MPAS_pool_get_array(velocitySolverPool, "elasticTimeStep", elasticTimeStep)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocityInitial", uVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocityInitial", vVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexU", airStressVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexV", airStressVertexV)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)
       call MPAS_pool_get_array(velocitySolverPool, "totalMassVertexfVertex", totalMassVertexfVertex)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)

       do iVertex = 1, nVerticesSolve

          if (solveVelocity(iVertex) == 1) then

             ! U equation
             leftMatrix(1,1) =  totalMassVertex(iVertex) / elasticTimeStep + oceanStressCoeff(iVertex) * cosOceanTurningAngle
             leftMatrix(1,2) = -totalMassVertexfVertex(iVertex) - &
                  oceanStressCoeff(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassVertexfVertex(iVertex))

             ! V equation
             leftMatrix(2,1) =  totalMassVertexfVertex(iVertex) + &
                  oceanStressCoeff(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassVertexfVertex(iVertex))
             leftMatrix(2,2) =  totalMassVertex(iVertex) / elasticTimeStep  + oceanStressCoeff(iVertex) * cosOceanTurningAngle

             ! right hand side of matrix solve
             rightHandSide(1) = stressDivergenceU(iVertex) + airStressVertexU(iVertex) + surfaceTiltForceU(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressU(iVertex) + (totalMassVertex(iVertex) * uVelocity(iVertex)) / elasticTimeStep

             rightHandSide(2) = stressDivergenceV(iVertex) + airStressVertexV(iVertex) + surfaceTiltForceV(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressV(iVertex) + (totalMassVertex(iVertex) * vVelocity(iVertex)) / elasticTimeStep

             ! solve the equation
             solutionDenominator = leftMatrix(1,1) * leftMatrix(2,2) - leftMatrix(1,2) * leftMatrix(2,1)

             uVelocity(iVertex) = (leftMatrix(2,2) * rightHandSide(1) - leftMatrix(1,2) * rightHandSide(2)) / solutionDenominator
             vVelocity(iVertex) = (leftMatrix(1,1) * rightHandSide(2) - leftMatrix(2,1) * rightHandSide(1)) / solutionDenominator

          endif ! solveVelocity

       enddo ! iVertex

       block => block % next
    end do

  end subroutine solve_velocity
  
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  solve_velocity_revised
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine solve_velocity_revised(domain)

    use cice_velocity_solver_constitutive_relation, only: &
         numericalInertiaCoefficient

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         icestatePool

    integer, dimension(:), pointer :: &
         solveVelocity

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         uVelocityInitial, &
         vVelocityInitial, &
         totalMassVertex, &
         totalMassVertexfVertex, &
         stressDivergenceU, &
         stressDivergenceV, &
         airStressVertexU, &
         airStressVertexV, &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         oceanStressU, &
         oceanStressV, &
         oceanStressCoeff

    real(kind=RKIND), pointer :: &
         dynamicsTimeStep

    real(kind=RKIND), dimension(2) :: &
         rightHandSide

    real(kind=RKIND), dimension(2,2) :: &
         leftMatrix

    real(kind=RKIND) :: &
         solutionDenominator

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))
 
       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocityInitial", uVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocityInitial", vVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexU", airStressVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexV", airStressVertexV)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)
       call MPAS_pool_get_array(velocitySolverPool, "totalMassVertexfVertex", totalMassVertexfVertex)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)
       call MPAS_pool_get_array(velocitySolverPool, "dynamicsTimeStep", dynamicsTimeStep)

       do iVertex = 1, nVerticesSolve

          if (solveVelocity(iVertex) == 1) then

             ! U equation
             leftMatrix(1,1) =  (numericalInertiaCoefficient + 1.0_RKIND) * (totalMassVertex(iVertex) / dynamicsTimeStep) &
                             +  oceanStressCoeff(iVertex) * cosOceanTurningAngle
             leftMatrix(1,2) = -totalMassVertexfVertex(iVertex) &
                             -  oceanStressCoeff(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassVertexfVertex(iVertex))

             ! V equation
             leftMatrix(2,1) =  totalMassVertexfVertex(iVertex) &
                             +  oceanStressCoeff(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassVertexfVertex(iVertex))
             leftMatrix(2,2) =  (numericalInertiaCoefficient + 1.0_RKIND) * (totalMassVertex(iVertex) / dynamicsTimeStep) &
                             +  oceanStressCoeff(iVertex) * cosOceanTurningAngle

             ! right hand side of matrix solve
             rightHandSide(1) = stressDivergenceU(iVertex) + airStressVertexU(iVertex) + surfaceTiltForceU(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressU(iVertex) + &
                  (totalMassVertex(iVertex) * (numericalInertiaCoefficient * uVelocity(iVertex) + uVelocityInitial(iVertex))) / dynamicsTimeStep

             rightHandSide(2) = stressDivergenceV(iVertex) + airStressVertexV(iVertex) + surfaceTiltForceV(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressV(iVertex) + &
                  (totalMassVertex(iVertex) * (numericalInertiaCoefficient * vVelocity(iVertex) + vVelocityInitial(iVertex))) / dynamicsTimeStep

             ! solve the equation
             solutionDenominator = leftMatrix(1,1) * leftMatrix(2,2) - leftMatrix(1,2) * leftMatrix(2,1)

             uVelocity(iVertex) = (leftMatrix(2,2) * rightHandSide(1) - leftMatrix(1,2) * rightHandSide(2)) / solutionDenominator
             vVelocity(iVertex) = (leftMatrix(1,1) * rightHandSide(2) - leftMatrix(2,1) * rightHandSide(1)) / solutionDenominator

          endif ! solveVelocity

       enddo ! iVertex

       block => block % next
    end do

  end subroutine solve_velocity_revised

!-----------------------------------------------------------------------
! Post sub-cycle
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  velocity_solver_post_subcycle
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date January 13th 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine velocity_solver_post_subcycle(domain)

    type(domain_type), intent(inout) :: &
         domain

    ! calculate the final divergence and shear for ridging
    call final_divergence_shear(domain)

    ! calculate principal stresses
    call principal_stresses_driver(domain)

    ! calculate final stress after subcycling
    call ocean_stress_final(domain)

  end subroutine velocity_solver_post_subcycle

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  final_divergence_shear
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date July 9th 2015
!> \details
!>  
!
!-----------------------------------------------------------------------
  
  subroutine final_divergence_shear(domain)

    use cice_velocity_solver_variational, only: &
         cice_final_divergence_shear_variational

    use cice_velocity_solver_weak, only: &
         cice_final_divergence_shear_weak

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block
    
    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)

       if (trim(config_stress_divergence_scheme) == "weak") then

          call cice_final_divergence_shear_weak(block)

       else if (trim(config_stress_divergence_scheme) == "variational") then

          call cice_final_divergence_shear_variational(block)

       endif

       block => block % next
    enddo

  end subroutine final_divergence_shear

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  principal_stresses_driver
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine principal_stresses_driver(domain)!{{{      

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         velocityWeakPool, &
         velocityVariationalPool

    real(kind=RKIND), dimension(:), pointer :: &
         replacementPressure, &
         principalStress1, &
         principalStress2

    real(kind=RKIND), dimension(:), pointer :: &
         stress11Weak, &
         stress22Weak, &
         stress12Weak

    real(kind=RKIND), dimension(:,:), pointer :: &
         stress11Var, &
         stress22Var, &
         stress12Var

    integer, pointer :: &
         nCellsSolve

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocityWeakPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocityVariationalPool)

       call MPAS_pool_get_array(velocitySolverPool, "replacementPressure", replacementPressure)
       call MPAS_pool_get_array(velocitySolverPool, "principalStress1", principalStress1)
       call MPAS_pool_get_array(velocitySolverPool, "principalStress2", principalStress2)

       ! calculate the principal stresses
       if (trim(config_stress_divergence_scheme) == "weak") then

          call MPAS_pool_get_array(velocityWeakPool, "stress11", stress11Weak)
          call MPAS_pool_get_array(velocityWeakPool, "stress22", stress22Weak)
          call MPAS_pool_get_array(velocityWeakPool, "stress12", stress12Weak)

          call principal_stresses(&
               nCellsSolve, &
               principalStress1, &
               principalStress2, &
               stress11Weak, &
               stress22Weak, &
               stress12Weak, &
               replacementPressure)

       else if (trim(config_stress_divergence_scheme) == "variational") then

          call MPAS_pool_get_array(velocityVariationalPool, "stress11", stress11Var)
          call MPAS_pool_get_array(velocityVariationalPool, "stress22", stress22Var)
          call MPAS_pool_get_array(velocityVariationalPool, "stress12", stress12Var)

          call principal_stresses(&
               nCellsSolve, &
               principalStress1, &
               principalStress2, &
               stress11Var(1,:), &
               stress22Var(1,:), &
               stress12Var(1,:), &
               replacementPressure)

       endif

       block => block % next
    enddo

  end subroutine principal_stresses_driver!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  principal_stresses
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine principal_stresses(&
       nPoints, &
       principalStress1, &
       principalStress2, &
       stress11, &
       stress22, &
       stress12, &
       replacementPressure)

    use cice_constants, only: &
         cicePuny

    integer, intent(in) :: &
         nPoints

    real(kind=RKIND), dimension(:), intent(out) :: &
         principalStress1, & !< Output: 
         principalStress2    !< Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         stress11, &         !< Input:
         stress22, &         !< Input:
         stress12, &         !< Input:
         replacementPressure !< Input:

    real(kind=RKIND) :: &
         sqrtContents

    integer :: &
         iPoint

    do iPoint = 1, nPoints

       if (replacementPressure(iPoint) > cicePuny) then

          sqrtContents = (stress11(iPoint) + stress22(iPoint))**2 - &
                         4.0_RKIND * stress11(iPoint) * stress22(iPoint) + &
                         4.0_RKIND * stress12(iPoint)**2

          principalStress1(iPoint) = 0.5_RKIND * (stress11(iPoint) + stress22(iPoint)) + 0.5_RKIND * sqrt(sqrtContents)
          principalStress2(iPoint) = 0.5_RKIND * (stress11(iPoint) + stress22(iPoint)) - 0.5_RKIND * sqrt(sqrtContents)

          principalStress1(iPoint) = principalStress1(iPoint) / replacementPressure(iPoint)
          principalStress2(iPoint) = principalStress2(iPoint) / replacementPressure(iPoint)

       else
          
          principalStress1(iPoint) = 1.0e30_RKIND
          principalStress2(iPoint) = 1.0e30_RKIND

       endif

    enddo ! iPoint

  end subroutine principal_stresses

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocean_stress_final
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine ocean_stress_final(domain)

    use cice_mesh, only: &
         cice_interpolate_vertex_to_cell

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         boundaryPool, &
         icestatePool, &
         velocitySolverPool

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaVertex, &
         fVertex, &
         uVelocity, &
         vVelocity, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         oceanStressU, &
         oceanStressV, &
         oceanStressCellU, &
         oceanStressCellV, &
         oceanStressCoeff

    integer, dimension(:), pointer :: &
         solveVelocity

    type(field1DReal), pointer :: &
         oceanStressUField, &
         oceanStressVField

    logical, pointer :: &
         configUseOceanStress

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    ! get ocean stress coefficient
    call ocean_stress_coefficient(domain)

    ! use stress config
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_ocean_stress", configUseOceanStress)

    ! get ocean stress on vertices
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellU", oceanStressCellU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellV", oceanStressCellV)

       if (configUseOceanStress) then

          call MPAS_pool_get_dimension(meshPool, "nVerticesSolve", nVerticesSolve)

          call MPAS_pool_get_array(meshPool, "fVertex", fVertex)

          call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)

          call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
          call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)
          call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)

          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) == 1) then

                oceanStressU(iVertex) = oceanStressCoeff(iVertex) * ((uOceanVelocityVertex(iVertex) - uVelocity(iVertex)) * cosOceanTurningAngle &
                                                                   - (vOceanVelocityVertex(iVertex) - vVelocity(iVertex)) * sinOceanTurningAngle * sign(1.0_RKIND,fVertex(iVertex)))
                oceanStressV(iVertex) = oceanStressCoeff(iVertex) * ((vOceanVelocityVertex(iVertex) - vVelocity(iVertex)) * cosOceanTurningAngle &
                                                                   + (uOceanVelocityVertex(iVertex) - uVelocity(iVertex)) * sinOceanTurningAngle * sign(1.0_RKIND,fVertex(iVertex)))

                oceanStressU(iVertex) = oceanStressU(iVertex) / iceAreaVertex(iVertex)
                oceanStressV(iVertex) = oceanStressV(iVertex) / iceAreaVertex(iVertex)

             else

                oceanStressU(iVertex) = 0.0_RKIND
                oceanStressV(iVertex) = 0.0_RKIND

             endif ! solveVelocity

          enddo ! iVertex

       else

          ! no ocean stress
          oceanStressU = 0.0_RKIND
          oceanStressV = 0.0_RKIND

       endif

       block => block % next
    end do

    ! get ocean stress on cells
    if (configUseOceanStress) then

       ! halo exchange ocean stress
       call mpas_timer_start("Halo")

       call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_field(velocitySolverPool, "oceanStressU", oceanStressUField)
       call MPAS_pool_get_field(velocitySolverPool, "oceanStressV", oceanStressVField)

       call MPAS_dmpar_exch_halo_field(oceanStressUField)
       call MPAS_dmpar_exch_halo_field(oceanStressVField)

       call mpas_timer_stop("Halo")

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

          call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
          call MPAS_pool_get_subpool(block % structs, "boundary", boundaryPool)
          call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
          call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

          call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)

          call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellU", oceanStressCellU)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellV", oceanStressCellV)

          ! interpolate ocean stress to cells
          call cice_interpolate_vertex_to_cell(meshPool, boundaryPool, oceanStressCellU, oceanStressU)
          call cice_interpolate_vertex_to_cell(meshPool, boundaryPool, oceanStressCellV, oceanStressV)

          ! multiply ocean stress back by ice area
          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) == 1) then

                oceanStressU(iVertex) = oceanStressU(iVertex) * iceAreaVertex(iVertex)
                oceanStressV(iVertex) = oceanStressV(iVertex) * iceAreaVertex(iVertex)

             endif ! solveVelocity

          enddo ! iVertex

          block => block % next
       end do

    else

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellU", oceanStressCellU)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellV", oceanStressCellV)

          oceanStressCellU = 0.0_RKIND
          oceanStressCellV = 0.0_RKIND

          block => block % next
       end do

    endif ! configUseOceanStress

  end subroutine ocean_stress_final

  !-----------------------------------------------------------------------

end module cice_velocity_solver
