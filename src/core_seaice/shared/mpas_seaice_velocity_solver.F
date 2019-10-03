!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_velocity_solver
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_velocity_solver

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_dmpar
  use mpas_timer
  use mpas_log, only: mpas_log_write

  implicit none

  private
  save

  public :: &
       seaice_init_velocity_solver, &
       seaice_run_velocity_solver

  ! velocity solver constants
  real(kind=RKIND), parameter, private :: &
       sinOceanTurningAngle = 0.0_RKIND, & ! northern hemisphere
       cosOceanTurningAngle = 1.0_RKIND, & ! northern hemisphere
       seaiceAreaMinimum = 0.001_RKIND, &
       seaiceMassMinimum = 0.01_RKIND

contains

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_velocity_solver
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_velocity_solver(&
       domain)!{{{

    use seaice_velocity_solver_constitutive_relation, only: &
         seaice_init_evp

    use seaice_velocity_solver_weak, only: &
         seaice_init_velocity_solver_weak

    use seaice_velocity_solver_variational, only: &
         seaice_init_velocity_solver_variational

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
         config_include_metric_terms, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch

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

    integer :: &
         ierr

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
       call seaice_init_evp(domain)

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

             call seaice_init_velocity_solver_weak(&
                  mesh, &
                  boundary, &
                  velocity_weak, &
                  config_rotate_cartesian_grid)

          else if (trim(config_stress_divergence_scheme) == "variational") then

             call seaice_init_velocity_solver_variational(&
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

             call MPAS_log_write("config_stress_divergence_scheme unknown: "//trim(config_stress_divergence_scheme), MPAS_LOG_CRIT)

          endif

          block => block % next
       enddo

    endif

    ! initialize the land ice shelve mask
    call init_ice_shelve_vertex_mask(domain)

    ! prep for aggregated halo exchanges
    call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
    if (config_aggregate_halo_exch) then

       ! create the velocity aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'velocityHaloExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create velocityHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'velocityHaloExchangeGroup', 'uVelocity', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add uVelocity to velocityHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'velocityHaloExchangeGroup', 'vVelocity', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add vVelocity to velocityHaloExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! create the iceAreaTotalMassExchangeGroup aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'iceAreaTotalMassExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create iceAreaTotalMassExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'iceAreaTotalMassExchangeGroup', 'iceAreaCellInitial', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add iceAreaCellInitial to iceAreaTotalMassExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'iceAreaTotalMassExchangeGroup', 'totalMassCell', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add totalMassCell to iceAreaTotalMassExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! create the solveVelocityExchangeGroup aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'solveVelocityExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create solveVelocityExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'solveVelocityExchangeGroup', 'solveVelocity', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add solveVelocity to solveVelocityExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! create the icePressureExchangeGroup aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'icePressureExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create icePressureExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'icePressureExchangeGroup', 'icePressure', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add icePressure to icePressureExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! create the airStressHaloExchangeGroup aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'airStressHaloExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create airStressHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'airStressHaloExchangeGroup', 'airStressCellU', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add airStressCellU to airStressHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'airStressHaloExchangeGroup', 'airStressCellV', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add airStressCellV to airStressHaloExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! create the seaSurfaceTiltHaloExchangeGroup aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'seaSurfaceTiltHaloExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create seaSurfaceTiltHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'seaSurfaceTiltHaloExchangeGroup', 'seaSurfaceTiltU', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add seaSurfaceTiltU to seaSurfaceTiltHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'seaSurfaceTiltHaloExchangeGroup', 'seaSurfaceTiltV', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add seaSurfaceTiltV to seaSurfaceTiltHaloExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! create the oceanStressHaloExchangeGroup aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'oceanStressHaloExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create oceanStressHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'oceanStressHaloExchangeGroup', 'oceanStressU', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add oceanStressU to oceanStressHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'oceanStressHaloExchangeGroup', 'oceanStressV', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add oceanStressV to oceanStressHaloExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! build reusable buffers
       call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
       if (config_reuse_halo_exch) then

          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "velocityHaloExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for velocityHaloExchangeGroup", MPAS_LOG_CRIT)
          endif
          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "iceAreaTotalMassExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for iceAreaTotalMassExchangeGroup", MPAS_LOG_CRIT)
          endif
          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "solveVelocityExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for solveVelocityExchangeGroup", MPAS_LOG_CRIT)
          endif
          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "icePressureExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for icePressureExchangeGroup", MPAS_LOG_CRIT)
          endif
          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "airStressHaloExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for airStressHaloExchangeGroup", MPAS_LOG_CRIT)
          endif
          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "seaSurfaceTiltHaloExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for seaSurfaceTiltHaloExchangeGroup", MPAS_LOG_CRIT)
          endif
          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "oceanStressHaloExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for oceanStressHaloExchangeGroup", MPAS_LOG_CRIT)
          endif

       endif

    endif ! config_aggregate_halo_exch

  end subroutine seaice_init_velocity_solver!}}}

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
!  seaice_run_velocity_solver
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_run_velocity_solver(domain, clock)!{{{

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

  end subroutine seaice_run_velocity_solver!}}}

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

    use seaice_constants, only: &
         seaiceDensityIce, &
         seaiceDensitySnow

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

          totalMassCell(iCell)  = iceVolumeCell(iCell)  * seaiceDensityIce + &
                                  snowVolumeCell(iCell) * seaiceDensitySnow

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

    use seaice_mesh, only: &
         seaice_interpolate_cell_to_vertex

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

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
         config_use_column_vertical_thermodynamics, &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch

    integer :: &
         ierr

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
    call seaice_load_balance_timers(domain, "vel prep before")

    call mpas_timer_start("ice area halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'iceAreaTotalMassExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for iceAreaTotalMassExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'iceAreaTotalMassExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for iceAreaTotalMassExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'iceAreaCellInitial')
          call MPAS_dmpar_field_halo_exch(domain, 'totalMassCell')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("ice area halo")

    call seaice_load_balance_timers(domain, "vel prep after")

    ! interpolate area and mass from cells to vertices
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

       call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)
       call MPAS_pool_get_array(icestatePool, "totalMassCell", totalMassCell)
       call MPAS_pool_get_array(icestatePool, "iceAreaCellInitial", iceAreaCellInitial)
       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

       call seaice_interpolate_cell_to_vertex(&
            meshPool, &
            iceAreaVertex, &
            iceAreaCellInitial)

       call seaice_interpolate_cell_to_vertex(&
            meshPool, &
            totalMassVertex, &
            totalMassCell)

       block => block % next
    end do

    ! calculate computational masks
    call stress_calculation_mask(domain)

    call velocity_calculation_mask(domain)

    ! halo exchange velocity mask
    call seaice_load_balance_timers(domain, "vel prep before")

    call mpas_timer_start("velocity mask halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'solveVelocityExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for solveVelocityExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'solveVelocityExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for solveVelocityExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'solveVelocity')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("velocity mask halo")

    call seaice_load_balance_timers(domain, "vel prep after")

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

          if (iceAreaCellInitial(iCell) > seaiceAreaMinimum .and. &
              totalMassCell(iCell) > seaiceMassMinimum .and. &
              landIceMask(iCell) == 0) then

             ! this cell has sufficient ice
             solveStress(iCell) = 1

          else

             ! test neighbouring cells to see if have sufficient ice
             do iCellOnCell = 1, nEdgesOnCell(iCell)

                iCellNeighbour = cellsOnCell(iCellOnCell,iCell)

                if (iceAreaCellInitial(iCellNeighbour) > seaiceAreaMinimum .and. &
                    totalMassCell(iCellNeighbour) > seaiceMassMinimum .and. &
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
              iceAreaVertex(iVertex) > seaiceAreaMinimum .and. &
              totalMassVertex(iVertex) > seaiceMassMinimum) then

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

    use seaice_mesh, only: &
         seaice_interpolate_cell_to_vertex

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

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
         iVertex, &
         ierr

    logical, pointer :: &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch

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
       call seaice_interpolate_cell_to_vertex(&
            meshPool, &
            uOceanVelocityVertex, &
            uOceanVelocity)

       call seaice_interpolate_cell_to_vertex(&
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
    call seaice_load_balance_timers(domain, "vel prep before")

    call mpas_timer_start("velocity halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'velocityHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for velocityHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'velocityHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for velocityHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'uVelocity')
          call MPAS_dmpar_field_halo_exch(domain, 'vVelocity')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("velocity halo")

    call seaice_load_balance_timers(domain, "vel prep after")

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

    use seaice_constants, only: &
         seaiceIceStrengthConstantHiblerP, &
         seaiceIceStrengthConstantHiblerC

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

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
         config_use_column_vertical_thermodynamics, &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    integer :: &
         iCell, &
         ierr

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

                icePressure(iCell) = seaiceIceStrengthConstantHiblerP * iceVolumeCell(iCell) * &
                     exp(-seaiceIceStrengthConstantHiblerC*(1.0_RKIND-iceAreaCell(iCell)))

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
    call seaice_load_balance_timers(domain, "vel prep before")

    call mpas_timer_start("ice strength halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'icePressureExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for icePressureExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'icePressureExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for icePressureExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'icePressure')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("ice strength halo")

    call seaice_load_balance_timers(domain, "vel prep after")

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

    use seaice_mesh, only: &
         seaice_interpolate_cell_to_vertex

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_vertical_thermodynamics, &
         config_use_air_stress, &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocitySolverPool

    real(kind=RKIND), dimension(:), pointer :: &
         airStressVertexU, &
         airStressVertexV, &
         airStressCellU, &
         airStressCellV

    integer :: &
         ierr

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
    call seaice_load_balance_timers(domain, "vel prep before")

    call mpas_timer_start("air stress halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'airStressHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for airStressHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'airStressHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for airStressHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'airStressCellU')
          call MPAS_dmpar_field_halo_exch(domain, 'airStressCellV')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("air stress halo")

    call seaice_load_balance_timers(domain, "vel prep after")

    ! interpolate air stress
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressCellV", airStressCellV)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexU", airStressVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexV", airStressVertexV)

       call seaice_interpolate_cell_to_vertex(&
            meshPool, &
            airStressVertexU, &
            airStressCellU)

       call seaice_interpolate_cell_to_vertex(&
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

    use seaice_mesh, only: &
         seaice_interpolate_cell_to_vertex

    use seaice_constants, only: &
         seaiceGravity

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

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
         iVertex, &
         ierr

    logical, pointer :: &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch

    ! halo exchange surface tilt
    call seaice_load_balance_timers(domain, "vel prep before")

    call mpas_timer_start("sea surface tilt halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'seaSurfaceTiltHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for seaSurfaceTiltHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'seaSurfaceTiltHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for seaSurfaceTiltHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'seaSurfaceTiltU')
          call MPAS_dmpar_field_halo_exch(domain, 'seaSurfaceTiltV')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("sea surface tilt halo")

    call seaice_load_balance_timers(domain, "vel prep after")

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
       call seaice_interpolate_cell_to_vertex(&
            meshPool, &
            seaSurfaceTiltVertexU, &
            seaSurfaceTiltU)

       call seaice_interpolate_cell_to_vertex(&
            meshPool, &
            seaSurfaceTiltVertexV, &
            seaSurfaceTiltV)

       ! calculate surface tilt from sea surface gradient
       do iVertex = 1, nVerticesSolve

          if (solveVelocity(iVertex) == 1) then

             surfaceTiltForceU(iVertex) = -seaiceGravity * totalMassVertex(iVertex) * seaSurfaceTiltVertexU(iVertex)
             surfaceTiltForceV(iVertex) = -seaiceGravity * totalMassVertex(iVertex) * seaSurfaceTiltVertexV(iVertex)

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

    use seaice_velocity_solver_weak, only: &
         seaice_internal_stress_weak

    use seaice_velocity_solver_variational, only: &
         seaice_internal_stress_variational

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, intent(in) :: &
         iElasticSubcycle !< Input: !! testing

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    logical, pointer :: &
         config_revised_evp, &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch

    integer :: &
         ierr

    call MPAS_pool_get_config(domain % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)
    call MPAS_pool_get_config(domain % configs, "config_revised_evp", config_revised_evp)

    ! calculate internal stresses
    if (trim(config_stress_divergence_scheme) == "weak") then

       call mpas_timer_start("Velocity solver internal stress")
       call seaice_internal_stress_weak(domain)
       call mpas_timer_stop("Velocity solver internal stress")

    else if (trim(config_stress_divergence_scheme) == "variational") then

       call mpas_timer_start("Velocity solver internal stress")
       call seaice_internal_stress_variational(domain)
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
    call seaice_load_balance_timers(domain, "vel before")

    call mpas_timer_start("Velocity solver halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'velocityHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for velocityHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'velocityHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for velocityHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'uVelocity')
          call MPAS_dmpar_field_halo_exch(domain, 'vVelocity')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("Velocity solver halo")

    call seaice_load_balance_timers(domain, "vel after")

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

    use seaice_constants, only: &
         seaiceDensitySeaWater, &
         seaiceIceOceanDragCoefficient

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

                oceanStressCoeff(iVertex) = seaiceIceOceanDragCoefficient * seaiceDensitySeaWater * iceAreaVertex(iVertex) * &
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

    use seaice_velocity_solver_constitutive_relation, only: &
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

    use seaice_velocity_solver_variational, only: &
         seaice_final_divergence_shear_variational

    use seaice_velocity_solver_weak, only: &
         seaice_final_divergence_shear_weak

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

          call seaice_final_divergence_shear_weak(block)

       else if (trim(config_stress_divergence_scheme) == "variational") then

          call seaice_final_divergence_shear_variational(block)

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

    use seaice_constants, only: &
         seaicePuny

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

    if (replacementPressure > seaicePuny) then

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

    use seaice_mesh, only: &
         seaice_interpolate_vertex_to_cell

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

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
         configUseOceanStress, &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex, &
         ierr

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
       call seaice_load_balance_timers(domain, "vel prep before")

       call mpas_timer_start("ocean stress halo")

       call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
       if (config_use_halo_exch) then

          call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
          if (config_aggregate_halo_exch) then

             ! aggregated halo exchange
             call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
             if (.not. config_reuse_halo_exch) then

                ! without reuse
                call mpas_dmpar_exch_group_full_halo_exch(domain, 'oceanStressHaloExchangeGroup', iErr=ierr)
                if (ierr /= MPAS_DMPAR_NOERR) then
                   call MPAS_log_write("failure to perform halo exchange for oceanStressHaloExchangeGroup", MPAS_LOG_CRIT)
                endif

             else

                ! with reuse
                call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'oceanStressHaloExchangeGroup', iErr=ierr)
                if (ierr /= MPAS_DMPAR_NOERR) then
                   call MPAS_log_write("failure to perform reuse halo exchange for oceanStressHaloExchangeGroup", MPAS_LOG_CRIT)
                endif

             endif ! config_reuse_halo_exch

          else

             ! no aggregated halo exchange
             call MPAS_dmpar_field_halo_exch(domain, 'oceanStressU')
             call MPAS_dmpar_field_halo_exch(domain, 'oceanStressV')

          endif ! config_aggregate_halo_exch

       endif ! config_use_halo_exch

       call mpas_timer_stop("ocean stress halo")

       call seaice_load_balance_timers(domain, "vel prep after")

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
          call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)

          ! interpolate ocean stress to cells
          call seaice_interpolate_vertex_to_cell(meshPool, boundaryPool, oceanStressCellU, oceanStressU)
          call seaice_interpolate_vertex_to_cell(meshPool, boundaryPool, oceanStressCellV, oceanStressV)

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

end module seaice_velocity_solver
