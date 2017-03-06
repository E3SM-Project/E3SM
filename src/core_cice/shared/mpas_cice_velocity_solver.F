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
         config_variational_basis, &
         config_wachspress_integration_type

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
         config_elastic_subcycle_number, &
         config_wachspress_integration_order

    ! set up the dynamically locked cell mask
    call dynamically_locked_cell_mask(domain)

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
          call MPAS_pool_get_config(block % configs, "config_wachspress_integration_type", config_wachspress_integration_type)
          call MPAS_pool_get_config(block % configs, "config_wachspress_integration_order", config_wachspress_integration_order)

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
                  config_variational_basis, &
                  config_wachspress_integration_type, &
                  config_wachspress_integration_order)

          else if (trim(config_stress_divergence_scheme) == "none") then

             continue

          else

             write(stderrUnit,*) "Error: config_stress_divergence_scheme unknown: ", trim(config_stress_divergence_scheme)
             call MPAS_dmpar_global_abort('MPAS-seaice: config_stress_divergence_scheme unknown')

          endif

          block => block % next
       enddo

    endif

    ! initialize the land ice shelve mask
    call init_ice_shelve_vertex_mask(domain)

  end subroutine cice_init_velocity_solver!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  dynamically_locked_cell_mask
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 6th April 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine dynamically_locked_cell_mask(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         meshPool, &
         boundaryPool

    integer, dimension(:), pointer :: &
         dynamicallyLockedCellsMask, &
         nEdgesOnCell, &
         interiorVertex

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    integer, pointer :: &
         nCells

    integer :: &
         iCell, &
         iVertexOnCell, &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "boundary", boundaryPool)

       call MPAS_pool_get_array(velocitySolverPool, "dynamicallyLockedCellsMask", dynamicallyLockedCellsMask)

       call MPAS_pool_get_array(boundaryPool, "interiorVertex", interiorVertex)

       call MPAS_pool_get_array(meshPool, "nEdgesOnCell", nEdgesOnCell)
       call MPAS_pool_get_array(meshPool, "verticesOnCell", verticesOnCell)

       call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

       do iCell = 1, nCells

          dynamicallyLockedCellsMask(iCell) = 1

          do iVertexOnCell = 1, nEdgesOnCell(iCell)

             iVertex = verticesOnCell(iVertexOnCell,iCell)

             if (interiorVertex(iVertex) == 1) then
                dynamicallyLockedCellsMask(iCell) = 0
                exit
             endif

          enddo ! iVertexOnCell

       enddo ! iCell

       block => block % next
    enddo

  end subroutine dynamically_locked_cell_mask

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_shelve_vertex_mask
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 9th April 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_shelve_vertex_mask(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         oceanCouplingPool, &
         meshPool

    integer, dimension(:), pointer :: &
         landIceMaskVertex, &
         landIceMask

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    integer, pointer :: &
         vertexDegree, &
         nVerticesSolve

    integer :: &
         iCell, &
         iVertex, &
         iCellOnVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(oceanCouplingPool, "landIceMask", landIceMask)
       call MPAS_pool_get_array(oceanCouplingPool, "landIceMaskVertex", landIceMaskVertex)

       call MPAS_pool_get_array(meshPool, "cellsOnVertex", cellsOnVertex)

       call MPAS_pool_get_dimension(block % dimensions, "vertexDegree", vertexDegree)
       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       landIceMaskVertex(:) = 0

       do iVertex = 1, nVerticesSolve

          do iCellOnVertex = 1, vertexDegree

             iCell = cellsOnVertex(iCellOnVertex,iVertex)

             if (landIceMask(iCell) == 1) then

                landIceMaskVertex(iVertex) = 1

             endif

          enddo ! iCellOnVertex

       enddo ! iVertex

       block => block % next
    enddo

  end subroutine init_ice_shelve_vertex_mask

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
    call mpas_timer_start("agg mass and area")
    call aggregate_mass_and_area(domain)
    call mpas_timer_stop("agg mass and area")

    ! calculate computational masks
    call mpas_timer_start("calc masks")
    call calculation_masks(domain)
    call mpas_timer_stop("calc masks")

    ! set new ice velocities to ocean velocity
    call mpas_timer_start("new ice vel")
    call new_ice_velocities(domain)
    call mpas_timer_stop("new ice vel")

    ! calculate the ice strength
    call mpas_timer_start("ice strength")
    call ice_strength(domain)
    call mpas_timer_stop("ice strength")

    ! calculate air stress
    call mpas_timer_start("air stress")
    call air_stress(domain)
    call mpas_timer_stop("air stress")

    ! calculate the coriolis force coefficient
    call mpas_timer_start("coliolis force coef")
    call coriolis_force_coefficient(domain)
    call mpas_timer_stop("coliolis force coef")

    ! calculate the ocean stress
    call mpas_timer_start("ocean stress")
    call ocean_stress(domain)
    call mpas_timer_stop("ocean stress")

    ! calculate the surface tilt force
    call mpas_timer_start("surface tilt")
    call surface_tilt(domain)
    call mpas_timer_stop("surface tilt")

    ! initialize subcycle variables
    call mpas_timer_start("init subcycle var")
    call init_subcycle_variables(domain)
    call mpas_timer_stop("init subcycle var")

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

    integer, pointer :: &
         nCellsSolve

    integer :: &
         iCell

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nCells", nCellsSolve)

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

       do iCell = 1, nCellsSolve

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
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaVertex, &
         iceAreaCellInitial, &
         iceAreaCell, &
         totalMassVertex, &
         totalMassCell

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_vertical_thermodynamics

    ! set initial ice area if not have column physics
    call MPAS_pool_get_config(domain % configs, "config_use_column_package", config_use_column_package)
    call MPAS_pool_get_config(domain % configs, "config_use_column_vertical_thermodynamics", &
                                                 config_use_column_vertical_thermodynamics)

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
    call mpas_timer_start("ice area halo")

    call MPAS_dmpar_field_halo_exch(domain, 'iceAreaCellInitial')
    call MPAS_dmpar_field_halo_exch(domain, 'totalMassCell')

    call mpas_timer_stop("ice area halo")

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
    call mpas_timer_start("velocity mask halo")

    call MPAS_dmpar_field_halo_exch(domain, 'solveVelocity')

    call mpas_timer_stop("velocity mask halo")

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
         icestatePool, &
         oceanCouplingPool

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
         nEdgesOnCell, &
         landIceMask

    integer, dimension(:,:), pointer :: &
         cellsOnCell

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)

       call MPAS_pool_get_array(meshPool, "nEdgesOnCell", nEdgesOnCell)
       call MPAS_pool_get_array(meshPool, "cellsOnCell", cellsOnCell)

       call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)

       call MPAS_pool_get_array(icestatePool, "iceAreaCellInitial", iceAreaCellInitial)
       call MPAS_pool_get_array(icestatePool, "totalMassCell", totalMassCell)

       call MPAS_pool_get_array(oceanCouplingPool, "landIceMask", landIceMask)

       do iCell = 1, nCells

          solveStress(iCell) = 0

          if (iceAreaCellInitial(iCell) > ciceAreaMinimum .and. &
              totalMassCell(iCell) > ciceMassMinimum .and. &
              landIceMask(iCell) == 0) then

             ! this cell has sufficient ice
             solveStress(iCell) = 1

          else

             ! test neighbouring cells to see if have sufficient ice
             do iCellOnCell = 1, nEdgesOnCell(iCell)

                iCellNeighbour = cellsOnCell(iCellOnCell,iCell)

                if (iceAreaCellInitial(iCellNeighbour) > ciceAreaMinimum .and. &
                    totalMassCell(iCellNeighbour) > ciceMassMinimum .and. &
                    landIceMask(iCellNeighbour) == 0) then

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
         boundaryPool, &
         oceanCouplingPool

    integer, dimension(:), pointer :: &
         solveVelocity

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaVertex, &
         totalMassVertex

    integer, dimension(:), pointer :: &
         interiorVertex, &
         landIceMaskVertex

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
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)

       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)

       call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)
       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

       call MPAS_pool_get_array(boundaryPool, "interiorVertex", interiorVertex)

       call MPAS_pool_get_array(oceanCouplingPool, "landIceMaskVertex", landIceMaskVertex)

       do iVertex = 1, nVerticesSolve

          solveVelocity(iVertex) = 0

          if (interiorVertex(iVertex) == 1 .and. &
              landIceMaskVertex(iVertex) == 0 .and. &
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
    call mpas_timer_start("velocity halo")

    call MPAS_dmpar_field_halo_exch(domain, 'uVelocity')
    call MPAS_dmpar_field_halo_exch(domain, 'vVelocity')

    call mpas_timer_stop("velocity halo")

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
       call MPAS_pool_get_config(block % configs, "config_use_column_vertical_thermodynamics", &
                                                   config_use_column_vertical_thermodynamics)

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
    call mpas_timer_start("ice strength halo")

    call MPAS_dmpar_field_halo_exch(domain, 'icePressure')

    call mpas_timer_stop("ice strength halo")

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

    ! calculate the air stress
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_package", config_use_column_package)
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_vertical_thermodynamics", &
                                                             config_use_column_vertical_thermodynamics)

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
    call mpas_timer_start("air stress halo")

    call MPAS_dmpar_field_halo_exch(domain, 'airStressCellU')
    call MPAS_dmpar_field_halo_exch(domain, 'airStressCellV')

    call mpas_timer_stop("air stress halo")

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

    logical, pointer :: &
         configUseSurfaceTilt, &
         configGeostrophicSurfaceTilt

    call MPAS_pool_get_config(domain % configs, "config_use_surface_tilt", configUseSurfaceTilt)
    call MPAS_pool_get_config(domain % configs, "config_geostrophic_surface_tilt", configGeostrophicSurfaceTilt)

    if (configUseSurfaceTilt) then

       if (configGeostrophicSurfaceTilt) then

          call surface_tilt_geostrophic(domain)

       else

          call surface_tilt_ssh_gradient(domain)

       endif

    else

       call no_surface_tilt(domain)

    endif

  end subroutine surface_tilt

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  surface_tilt_geostrophic
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine surface_tilt_geostrophic(domain)

    use cice_constants, only: &
         ciceGravity

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

       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)

       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

       call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
       call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)
       call MPAS_pool_get_array(meshPool, "fVertex", fVertex)

       ! calculate surface tilt force from geostrophic currents
       do iVertex = 1, nVerticesSolve

          if (solveVelocity(iVertex) == 1) then

             surfaceTiltForceU(iVertex) = -fVertex(iVertex) * totalMassVertex(iVertex) * vOceanVelocityVertex(iVertex)
             surfaceTiltForceV(iVertex) =  fVertex(iVertex) * totalMassVertex(iVertex) * uOceanVelocityVertex(iVertex)

          else

             surfaceTiltForceU(iVertex) = 0.0_RKIND
             surfaceTiltForceV(iVertex) = 0.0_RKIND

          endif ! solveVelocity

       enddo ! iVertex

       block => block % next
    end do

  end subroutine surface_tilt_geostrophic

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  surface_tilt_ssh_gradient
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine surface_tilt_ssh_gradient(domain)

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    use cice_constants, only: &
         ciceGravity

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         icestatePool, &
         velocitySolverPool, &
         meshPool, &
         oceanCouplingPool

    real(kind=RKIND), dimension(:), pointer:: &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         totalMassVertex, &
         seaSurfaceTiltU, &
         seaSurfaceTiltV, &
         seaSurfaceTiltVertexU, &
         seaSurfaceTiltVertexV

    integer, dimension(:), pointer :: &
         solveVelocity

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    ! halo exchange surface tilt
    call mpas_timer_start("sea surface tilt halo")

    call MPAS_dmpar_field_halo_exch(domain, 'seaSurfaceTiltU')
    call MPAS_dmpar_field_halo_exch(domain, 'seaSurfaceTiltV')

    call mpas_timer_stop("sea surface tilt halo")

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)

       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

       call MPAS_pool_get_array(oceanCouplingPool, "seaSurfaceTiltU", seaSurfaceTiltU)
       call MPAS_pool_get_array(oceanCouplingPool, "seaSurfaceTiltV", seaSurfaceTiltV)

       call MPAS_pool_get_array(velocitySolverPool, "seaSurfaceTiltVertexU", seaSurfaceTiltVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "seaSurfaceTiltVertexV", seaSurfaceTiltVertexV)

       ! interpolate sea surface tilt from cells to vertices
       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            seaSurfaceTiltVertexU, &
            seaSurfaceTiltU)

       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            seaSurfaceTiltVertexV, &
            seaSurfaceTiltV)

       ! calculate surface tilt from sea surface gradient
       do iVertex = 1, nVerticesSolve

          if (solveVelocity(iVertex) == 1) then

             surfaceTiltForceU(iVertex) = -ciceGravity * totalMassVertex(iVertex) * seaSurfaceTiltVertexU(iVertex)
             surfaceTiltForceV(iVertex) = -ciceGravity * totalMassVertex(iVertex) * seaSurfaceTiltVertexV(iVertex)

          else

             surfaceTiltForceU(iVertex) = 0.0_RKIND
             surfaceTiltForceV(iVertex) = 0.0_RKIND

          endif ! solveVelocity

       enddo ! iVertex

       block => block % next
    end do

  end subroutine surface_tilt_ssh_gradient

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  no_surface_tilt
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine no_surface_tilt(domain)

    use cice_constants, only: &
         ciceGravity

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool

    real(kind=RKIND), dimension(:), pointer:: &
         surfaceTiltForceU, &
         surfaceTiltForceV

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)

       ! no surface tilt
       surfaceTiltForceU = 0.0_RKIND
       surfaceTiltForceV = 0.0_RKIND

       block => block % next
    end do

  end subroutine no_surface_tilt

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
         velocityVariationalPool, &
         velocityWeakPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11Var, &
         strain22Var, &
         strain12Var, &
         stress11Var, &
         stress22Var, &
         stress12Var

    real(kind=RKIND), dimension(:), pointer :: &
         strain11Weak, &
         strain22Weak, &
         strain12Weak, &
         stress11Weak, &
         stress22Weak, &
         stress12Weak

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

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    integer :: &
         iCell, &
         iVertex

    call MPAS_pool_get_config(domain % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

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

       if (trim(config_stress_divergence_scheme) == "variational") then
          ! variational

          call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocityVariationalPool)

          ! strains
          call MPAS_pool_get_array(velocityVariationalPool, "strain11", strain11Var)
          call MPAS_pool_get_array(velocityVariationalPool, "strain22", strain22Var)
          call MPAS_pool_get_array(velocityVariationalPool, "strain12", strain12Var)

          strain11Var = 0.0_RKIND
          strain22Var = 0.0_RKIND
          strain12Var = 0.0_RKIND

          ! stresses
          call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)

          call MPAS_pool_get_array(velocityVariationalPool, "stress11", stress11Var)
          call MPAS_pool_get_array(velocityVariationalPool, "stress22", stress22Var)
          call MPAS_pool_get_array(velocityVariationalPool, "stress12", stress12Var)

          call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

          do iCell = 1, nCells

             if (solveStress(iCell) /= 1) then

                stress11Var(:,iCell) = 0.0_RKIND
                stress22Var(:,iCell) = 0.0_RKIND
                stress12Var(:,iCell) = 0.0_RKIND

             endif

          enddo ! iCell

       else if (trim(config_stress_divergence_scheme) == "weak") then

          call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocityWeakPool)

          ! strains
          call MPAS_pool_get_array(velocityWeakPool, "strain11", strain11Weak)
          call MPAS_pool_get_array(velocityWeakPool, "strain22", strain22Weak)
          call MPAS_pool_get_array(velocityWeakPool, "strain12", strain12Weak)

          strain11Weak = 0.0_RKIND
          strain22Weak = 0.0_RKIND
          strain12Weak = 0.0_RKIND

          ! stresses
          call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)

          call MPAS_pool_get_array(velocityWeakPool, "stress11", stress11Weak)
          call MPAS_pool_get_array(velocityWeakPool, "stress22", stress22Weak)
          call MPAS_pool_get_array(velocityWeakPool, "stress12", stress12Weak)

          call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

          do iCell = 1, nCells

             if (solveStress(iCell) /= 1) then

                stress11Weak(iCell) = 0.0_RKIND
                stress22Weak(iCell) = 0.0_RKIND
                stress12Weak(iCell) = 0.0_RKIND

             endif

          enddo ! iCell

       endif ! config_stress_divergence_scheme

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
    call mpas_timer_start("ocn stress coef")
    call ocean_stress_coefficient(domain)
    call mpas_timer_stop("ocn stress coef")

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
    call mpas_timer_start("Velocity solver halo")

    call MPAS_dmpar_field_halo_exch(domain, 'uVelocity')
    call MPAS_dmpar_field_halo_exch(domain, 'vVelocity')

    call mpas_timer_stop("Velocity solver halo")

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

       !$omp parallel do default(shared) private(iVertex, leftMatrix, rightHandSide, solutionDenominator)
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
                  oceanStressCoeff(iVertex) * oceanStressU(iVertex) + &
                  (totalMassVertex(iVertex) * uVelocity(iVertex)) / elasticTimeStep

             rightHandSide(2) = stressDivergenceV(iVertex) + airStressVertexV(iVertex) + surfaceTiltForceV(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressV(iVertex) + &
                  (totalMassVertex(iVertex) * vVelocity(iVertex)) / elasticTimeStep

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
                  (totalMassVertex(iVertex) * (numericalInertiaCoefficient * uVelocity(iVertex) + &
                  uVelocityInitial(iVertex))) / dynamicsTimeStep

             rightHandSide(2) = stressDivergenceV(iVertex) + airStressVertexV(iVertex) + surfaceTiltForceV(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressV(iVertex) + &
                  (totalMassVertex(iVertex) * (numericalInertiaCoefficient * vVelocity(iVertex) + &
                  vVelocityInitial(iVertex))) / dynamicsTimeStep

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
    call mpas_timer_start("final div shear")
    call final_divergence_shear(domain)
    call mpas_timer_stop("final div shear")

    ! calculate principal stresses
    call mpas_timer_start("principal stress")
    call principal_stresses_driver(domain)
    call mpas_timer_stop("principal stress")

    ! calculate final stress after subcycling
    call mpas_timer_start("ocn stress final")
    call ocean_stress_final(domain)
    call mpas_timer_stop("ocn stress final")

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
         velocityWeakPool, &
         velocityVariationalPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         replacementPressureWeak, &
         principalStress1Weak, &
         principalStress2Weak

    real(kind=RKIND), dimension(:,:), pointer :: &
         replacementPressureVar, &
         principalStress1Var, &
         principalStress2Var

    real(kind=RKIND), dimension(:), pointer :: &
         stress11Weak, &
         stress22Weak, &
         stress12Weak

    real(kind=RKIND), dimension(:,:), pointer :: &
         stress11Var, &
         stress22Var, &
         stress12Var

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, pointer :: &
         nCellsSolve

    integer :: &
         iCell, &
         iVertexOnCell

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! calculate the principal stresses
       if (trim(config_stress_divergence_scheme) == "weak") then

          call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocityWeakPool)

          call MPAS_pool_get_array(velocityWeakPool, "stress11", stress11Weak)
          call MPAS_pool_get_array(velocityWeakPool, "stress22", stress22Weak)
          call MPAS_pool_get_array(velocityWeakPool, "stress12", stress12Weak)
          call MPAS_pool_get_array(velocityWeakPool, "replacementPressure", replacementPressureWeak)
          call MPAS_pool_get_array(velocityWeakPool, "principalStress1", principalStress1Weak)
          call MPAS_pool_get_array(velocityWeakPool, "principalStress2", principalStress2Weak)

          do iCell = 1, nCellsSolve

             call principal_stresses(&
                  principalStress1Weak(iCell), &
                  principalStress2Weak(iCell), &
                  stress11Weak(iCell), &
                  stress22Weak(iCell), &
                  stress12Weak(iCell), &
                  replacementPressureWeak(iCell))

          enddo ! iCell

       else if (trim(config_stress_divergence_scheme) == "variational") then

          call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocityVariationalPool)
          call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

          call MPAS_pool_get_array(meshPool, "nEdgesOnCell", nEdgesOnCell)

          call MPAS_pool_get_array(velocityVariationalPool, "stress11", stress11Var)
          call MPAS_pool_get_array(velocityVariationalPool, "stress22", stress22Var)
          call MPAS_pool_get_array(velocityVariationalPool, "stress12", stress12Var)
          call MPAS_pool_get_array(velocityVariationalPool, "replacementPressure", replacementPressureVar)
          call MPAS_pool_get_array(velocityVariationalPool, "principalStress1", principalStress1Var)
          call MPAS_pool_get_array(velocityVariationalPool, "principalStress2", principalStress2Var)

          do iCell = 1, nCellsSolve
             do iVertexOnCell = 1, nEdgesOnCell(iCell)

                call principal_stresses(&
                     principalStress1Var(iVertexOnCell,iCell), &
                     principalStress2Var(iVertexOnCell,iCell), &
                     stress11Var(iVertexOnCell,iCell), &
                     stress22Var(iVertexOnCell,iCell), &
                     stress12Var(iVertexOnCell,iCell), &
                     replacementPressureVar(iVertexOnCell,iCell))

             enddo ! iVertexOnCell
          enddo ! iCell

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
       principalStress1, &
       principalStress2, &
       stress11, &
       stress22, &
       stress12, &
       replacementPressure)

    use cice_constants, only: &
         cicePuny

    real(kind=RKIND), intent(out) :: &
         principalStress1, & !< Output:
         principalStress2    !< Output:

    real(kind=RKIND), intent(in) :: &
         stress11, &         !< Input:
         stress22, &         !< Input:
         stress12, &         !< Input:
         replacementPressure !< Input:

    real(kind=RKIND) :: &
         sqrtContents

    if (replacementPressure > cicePuny) then

       sqrtContents = (stress11 + stress22)**2 - &
                      4.0_RKIND * stress11 * stress22 + &
                      4.0_RKIND * stress12**2

       principalStress1 = 0.5_RKIND * (stress11 + stress22) + 0.5_RKIND * sqrt(sqrtContents)
       principalStress2 = 0.5_RKIND * (stress11 + stress22) - 0.5_RKIND * sqrt(sqrtContents)

       principalStress1 = principalStress1 / replacementPressure
       principalStress2 = principalStress2 / replacementPressure

    else

       principalStress1 = 1.0e30_RKIND
       principalStress2 = 1.0e30_RKIND

    endif

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

                oceanStressU(iVertex) = oceanStressCoeff(iVertex) * &
                     ((uOceanVelocityVertex(iVertex) - uVelocity(iVertex)) * cosOceanTurningAngle &
                    - (vOceanVelocityVertex(iVertex) - vVelocity(iVertex)) * sinOceanTurningAngle * &
                     sign(1.0_RKIND,fVertex(iVertex)))
                oceanStressV(iVertex) = oceanStressCoeff(iVertex) * &
                     ((vOceanVelocityVertex(iVertex) - vVelocity(iVertex)) * cosOceanTurningAngle &
                    + (uOceanVelocityVertex(iVertex) - uVelocity(iVertex)) * sinOceanTurningAngle * &
                     sign(1.0_RKIND,fVertex(iVertex)))

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
       call mpas_timer_start("ocean stress halo")

       call MPAS_dmpar_field_halo_exch(domain, 'oceanStressU')
       call MPAS_dmpar_field_halo_exch(domain, 'oceanStressV')

       call mpas_timer_stop("ocean stress halo")

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
