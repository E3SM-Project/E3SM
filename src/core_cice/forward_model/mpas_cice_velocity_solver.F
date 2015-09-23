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

  implicit none

  private
  save

  public :: &
       cice_init_velocity_solver, &
       cice_run_velocity_solver
  
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

    use cice_velocity_solver_shared, only: &
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
         config_dt

    ! set things even if we dont want velocity
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolver)
       call MPAS_pool_get_array(velocitySolver, "dynamicsTimeStep", dynamicsTimeStep)

       dynamicsTimeStep = config_dt

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

  subroutine cice_run_velocity_solver(&
       domain, &
       clock, &
       dt)!{{{

    use cice_debug

    type (domain_type), intent(inout) :: &
         domain !< Input/Output: 

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input: 

    real(kind=RKIND), intent(in) :: &
         dt !< Input: 

    logical, pointer :: &
         config_use_velocity_solver

    ! determine if velocity solver switched on
    call MPAS_pool_get_config(domain % configs, "config_use_velocity_solver", config_use_velocity_solver)

    if (config_use_velocity_solver) then

       call cice_debug_write_out_fields(domain, "cice_run_velocity_solver", 1)

       ! pre subcycle
       call velocity_solver_pre_subcycle(domain)

       ! subcycle the dynamics
       call subcycle_velocity_solver(domain, clock)

       ! post subcycle
       call velocity_solver_post_subcycle(domain)

       ! geographical output
       call geographical_velocity_diagnostics(domain)
       
       call cice_debug_write_out_fields(domain, "cice_run_velocity_solver", 2)

    endif ! config_use_velocity_solver

  end subroutine cice_run_velocity_solver!}}}

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

    use cice_velocity_solver_shared, only: &
         cice_ice_strength, &
         cice_air_stress, &
         cice_coriolis_force_coefficient, &
         cice_ocean_stress, &
         cice_surface_tilt

    type(domain_type), intent(inout) :: &
         domain

    ! aggregate categories for area and volume into total mass
    call aggregate_mass_and_area(domain)

    ! calculate computational masks
    call calculation_masks(domain)

    ! set new ice velocities to ocean velocity
    call new_ice_velocities(domain)

    ! calculate the ice strength
    call cice_ice_strength(domain)

    ! calculate air stress
    call cice_air_stress(domain)

    ! calculate the coriolis force coefficient
    call cice_coriolis_force_coefficient(domain)

    ! calculate the ocean stress
    call cice_ocean_stress(domain)

    ! calculate the surface tilt force
    call cice_surface_tilt(domain)

  end subroutine velocity_solver_pre_subcycle

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

    use cice_velocity_solver_shared, only: &
         cice_ocean_stress_final

    type(domain_type), intent(inout) :: &
         domain

    ! calculate the final divergence and shear for ridging
    call final_divergence_shear(domain)

    ! calculate principal stresses
    call principal_stresses_driver(domain)

    ! calculate final stress after subcycling
    call cice_ocean_stress_final(domain)

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

    use cice_velocity_solver_shared, only: &
         cice_principal_stresses         

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

          call cice_principal_stresses(&
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

          call cice_principal_stresses(&
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

    use cice_velocity_solver_shared, only: &
         cice_ocean_stress_coefficient, &
         cice_solve_velocity, &
         cice_solve_velocity_revised

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

       call cice_internal_stress_weak(domain)

    else if (trim(config_stress_divergence_scheme) == "variational") then

       call cice_internal_stress_variational(domain)

    endif

    ! ocean stress coefficient
    call cice_ocean_stress_coefficient(domain)

    ! solve for velocity
    if (.not. config_revised_evp) then

       call cice_solve_velocity(domain)

    else

       call cice_solve_velocity_revised(domain)

    endif

    ! halo exchange
    call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolverPool)

    call MPAS_pool_get_field(velocitySolverPool, "uVelocity", uVelocity)
    call MPAS_pool_get_field(velocitySolverPool, "vVelocity", vVelocity)

    call MPAS_dmpar_exch_halo_field(uVelocity)
    call MPAS_dmpar_exch_halo_field(vVelocity)

  end subroutine single_subcycle_velocity_solver!}}}

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
    call MPAS_pool_get_subpool(domain % blocklist % structs, "tracers", tracersPool)

    call MPAS_pool_get_field(tracersPool, "iceAreaCategory", iceAreaCategoryField, 1)
    call MPAS_pool_get_field(tracersPool, "iceVolumeCategory", iceVolumeCategoryField, 1)
    call MPAS_pool_get_field(tracersPool, "snowVolumeCategory", snowVolumeCategoryField, 1)

    call MPAS_dmpar_exch_halo_field(iceAreaCategoryField)
    call MPAS_dmpar_exch_halo_field(iceVolumeCategoryField)
    call MPAS_dmpar_exch_halo_field(snowVolumeCategoryField)

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
         meshPool, &
         velocitySolverPool

    type(field1DReal), pointer :: &
         iceAreaCellInitialField

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaVertex, &
         iceAreaCellInitial, &
         totalMassVertex, &
         totalMassCell

    type(field1DInteger), pointer :: &
         solveVelocityField

    ! halo exchange of initial ice area of cell
    call MPAS_pool_get_subpool(domain % blocklist % structs, "icestate", icestatePool)
    call MPAS_pool_get_field(icestatePool, "iceAreaCellInitial", iceAreaCellInitialField)
    call MPAS_dmpar_exch_halo_field(iceAreaCellInitialField)

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
    call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolverPool)
    call MPAS_pool_get_field(velocitySolverPool, "solveVelocity", solveVelocityField)
    call MPAS_dmpar_exch_halo_field(solveVelocityField)

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

    use cice_velocity_solver_shared, only: &
         ciceAreaMinimum, &
         ciceMassMinimum

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

    use cice_velocity_solver_shared, only: &
         ciceAreaMinimum, &
         ciceMassMinimum

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
    call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolverPool)
    call MPAS_pool_get_field(velocitySolverPool, "uVelocity", uVelocityField)
    call MPAS_pool_get_field(velocitySolverPool, "vVelocity", vVelocityField)
    call MPAS_dmpar_exch_halo_field(uVelocityField)
    call MPAS_dmpar_exch_halo_field(vVelocityField)

  end subroutine new_ice_velocities!}}}
    
  !-----------------------------------------------------------------------
  ! velocity diagnostics
  !-----------------------------------------------------------------------

  !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  !
  !  geographical_velocity_diagnostics
  !
  !> \brief
  !> \author Adrian K. Turner, LANL
  !> \date
  !> \details
  !>
  !
  !-----------------------------------------------------------------------

  subroutine geographical_velocity_diagnostics(domain)

    use cice_mesh, only: &
         cice_latlon_vector_rotation_backward

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocitySolverPool, &
         velocityGeographicalPool

    ! mesh quantities
    real(kind=RKIND), dimension(:), pointer :: &
         latVertex, &
         lonVertex, &
         xVertex, &
         yVertex, &
         zVertex

    ! other quantities
    real(kind=RKIND), pointer :: &
         sphereRadius

    logical, pointer :: &
         rotateCartesianGrid

    ! non-geographical arrays
    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         stressDivergenceU, &
         stressDivergenceV, &
         airStressVertexU, &
         airStressVertexV, &
         oceanStressU, &
         oceanStressV, &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex

    ! geographical arrays
    real(kind=RKIND), dimension(:), pointer :: &
         uVelocityGeo, &
         vVelocityGeo, &
         stressDivergenceUGeo, &
         stressDivergenceVGeo, &
         airStressVertexUGeo, &
         airStressVertexVGeo, &
         oceanStressUGeo, &
         oceanStressVGeo, &
         surfaceTiltForceUGeo, &
         surfaceTiltForceVGeo, &
         uOceanVelocityVertexGeo, &
         vOceanVelocityVertexGeo

    integer, pointer :: &
         nVertices

    integer :: &
         iVertex

    ! get configs
    call MPAS_pool_get_config(domain % configs, "config_rotate_cartesian_grid", rotateCartesianGrid)

    block => domain % blocklist
    do while (associated(block))

       ! get pools
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_geographical", velocityGeographicalPool)

       ! mesh arrays
       call MPAS_pool_get_config(meshPool, "sphere_radius", sphereRadius)
       call MPAS_pool_get_array(meshPool, "latVertex", latVertex)
       call MPAS_pool_get_array(meshPool, "lonVertex", lonVertex)
       call MPAS_pool_get_array(meshPool, "xVertex", xVertex)
       call MPAS_pool_get_array(meshPool, "yVertex", yVertex)
       call MPAS_pool_get_array(meshPool, "zVertex", zVertex)

       call MPAS_pool_get_dimension(meshPool, "nVertices", nVertices)

       ! get non-geographical arrays
       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexU", airStressVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexV", airStressVertexV)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)
       call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
       call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)

       ! get geographical arrays
       call MPAS_pool_get_array(velocityGeographicalPool, "uVelocityGeo", uVelocityGeo)
       call MPAS_pool_get_array(velocityGeographicalPool, "vVelocityGeo", vVelocityGeo)
       call MPAS_pool_get_array(velocityGeographicalPool, "stressDivergenceUGeo", stressDivergenceUGeo)
       call MPAS_pool_get_array(velocityGeographicalPool, "stressDivergenceVGeo", stressDivergenceVGeo)
       call MPAS_pool_get_array(velocityGeographicalPool, "airStressVertexUGeo", airStressVertexUGeo)
       call MPAS_pool_get_array(velocityGeographicalPool, "airStressVertexVGeo", airStressVertexVGeo)
       call MPAS_pool_get_array(velocityGeographicalPool, "oceanStressUGeo", oceanStressUGeo)
       call MPAS_pool_get_array(velocityGeographicalPool, "oceanStressVGeo", oceanStressVGeo)
       call MPAS_pool_get_array(velocityGeographicalPool, "surfaceTiltForceUGeo", surfaceTiltForceUGeo)
       call MPAS_pool_get_array(velocityGeographicalPool, "surfaceTiltForceVGeo", surfaceTiltForceVGeo)
       call MPAS_pool_get_array(velocityGeographicalPool, "uOceanVelocityVertexGeo", uOceanVelocityVertexGeo)
       call MPAS_pool_get_array(velocityGeographicalPool, "vOceanVelocityVertexGeo", vOceanVelocityVertexGeo)

       ! perform rotations
       do iVertex = 1, nVertices

          ! ice velocity
          call cice_latlon_vector_rotation_backward(&
               uVelocityGeo(iVertex), &
               vVelocityGeo(iVertex), &
               uVelocity(iVertex), &
               vVelocity(iVertex), &
               latVertex(iVertex), &
               lonVertex(iVertex), &
               xVertex(iVertex), &
               yVertex(iVertex), &
               zVertex(iVertex), &
               sphereRadius, &
               rotateCartesianGrid)

          ! stress divergence
          call cice_latlon_vector_rotation_backward(&
               stressDivergenceUGeo(iVertex), &
               stressDivergenceVGeo(iVertex), &
               stressDivergenceU(iVertex), &
               stressDivergenceV(iVertex), &
               latVertex(iVertex), &
               lonVertex(iVertex), &
               xVertex(iVertex), &
               yVertex(iVertex), &
               zVertex(iVertex), &
               sphereRadius, &
               rotateCartesianGrid)

          ! air stress
          call cice_latlon_vector_rotation_backward(&
               airStressVertexUGeo(iVertex), &
               airStressVertexVGeo(iVertex), &
               airStressVertexU(iVertex), &
               airStressVertexV(iVertex), &
               latVertex(iVertex), &
               lonVertex(iVertex), &
               xVertex(iVertex), &
               yVertex(iVertex), &
               zVertex(iVertex), &
               sphereRadius, &
               rotateCartesianGrid)

          ! ocean stress
          call cice_latlon_vector_rotation_backward(&
               oceanStressUGeo(iVertex), &
               oceanStressVGeo(iVertex), &
               oceanStressU(iVertex), &
               oceanStressV(iVertex), &
               latVertex(iVertex), &
               lonVertex(iVertex), &
               xVertex(iVertex), &
               yVertex(iVertex), &
               zVertex(iVertex), &
               sphereRadius, &
               rotateCartesianGrid)

          ! surface tilt
          call cice_latlon_vector_rotation_backward(&
               surfaceTiltForceUGeo(iVertex), &
               surfaceTiltForceVGeo(iVertex), &
               surfaceTiltForceU(iVertex), &
               surfaceTiltForceV(iVertex), &
               latVertex(iVertex), &
               lonVertex(iVertex), &
               xVertex(iVertex), &
               yVertex(iVertex), &
               zVertex(iVertex), &
               sphereRadius, &
               rotateCartesianGrid)

          ! ocean velocity
          call cice_latlon_vector_rotation_backward(&
               uOceanVelocityVertexGeo(iVertex), &
               vOceanVelocityVertexGeo(iVertex), &
               uOceanVelocityVertex(iVertex), &
               vOceanVelocityVertex(iVertex), &
               latVertex(iVertex), &
               lonVertex(iVertex), &
               xVertex(iVertex), &
               yVertex(iVertex), &
               zVertex(iVertex), &
               sphereRadius, &
               rotateCartesianGrid)

       enddo ! iVertex

       block => block % next
    enddo

  end subroutine geographical_velocity_diagnostics

  !-----------------------------------------------------------------------

end module cice_velocity_solver
