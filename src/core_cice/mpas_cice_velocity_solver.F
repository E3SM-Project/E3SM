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

    type(block_type), pointer :: &
         block

    type (MPAS_pool_type), pointer :: &
         tracers_fields

    type(field3DReal), pointer :: &
         iceAreaCategoryField, &
         iceVolumeCategoryField, &
         snowVolumeCategoryField

    logical, pointer :: &
         config_use_velocity_solver

    call MPAS_pool_get_config(domain % configs, "config_use_velocity_solver", config_use_velocity_solver)

    if (config_use_velocity_solver) then

       call cice_debug_write_out_fields(domain, "cice_run_velocity_solver", 1)

       ! halo updates
       call MPAS_pool_get_subpool(domain % blocklist % structs, "tracers", tracers_fields)

       call MPAS_pool_get_field(tracers_fields, "iceAreaCategory", iceAreaCategoryField, 1)
       call MPAS_pool_get_field(tracers_fields, "iceVolumeCategory", iceVolumeCategoryField, 1)
       call MPAS_pool_get_field(tracers_fields, "snowVolumeCategory", snowVolumeCategoryField, 1)

       call MPAS_dmpar_exch_halo_field(iceAreaCategoryField)
       call MPAS_dmpar_exch_halo_field(iceVolumeCategoryField)
       call MPAS_dmpar_exch_halo_field(snowVolumeCategoryField)

       ! pre subcycle
       block => domain % blocklist
       do while (associated(block))
          call velocity_solver_pre_subcycle(block)
          block => block % next
       end do

       ! subcycle the dynamics
       call subcycle_velocity_solver(&
            domain, &
            clock)

       ! post subcycle
       block => domain % blocklist
       do while (associated(block))
          call velocity_solver_post_subcycle(block)
          block => block % next
       end do

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

  subroutine velocity_solver_pre_subcycle(block)

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    use cice_velocity_solver_shared, only: &
         cice_ice_strength, &
         cice_air_stress, &
         cice_coriolis_force_coefficient, &
         cice_ocean_stress, &
         cice_surface_tilt

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         icestate, &
         tracers_aggregate, &
         velocity_solver, &
         boundary, &
         atmos_coupling, &
         ocean_coupling

    integer, pointer :: &
         nCells, &
         nVerticesSolve, &
         nCategories

    integer, dimension(:), pointer :: &
         solveStress, &
         solveVelocity, &
         solveVelocityPrevious, &
         interiorVertex

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         uVelocityInitial, &
         vVelocityInitial, & 
         iceAreaCell, &
         totalMassCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         iceAreaCellInitial, &
         openWaterArea, &
         iceAreaVertex, &
         totalMassVertex, &
         icePressure, &
         airStressCellU, &
         airStressCellV, &
         airStressVertexU, &
         airStressVertexV, &
         oceanStressU, &
         oceanStressV, &
         stressDivergenceU, &
         stressDivergenceV, &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         totalMassVertexfVertex, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         airDensity, &
         uAirVelocity, &
         vAirVelocity, &
         uOceanVelocity, &
         vOceanVelocity, &
         fVertex

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_vertical_thermodynamics

    call MPAS_pool_get_config(block % configs, "config_use_column_package", config_use_column_package)
    call MPAS_pool_get_config(block % configs, "config_use_column_vertical_thermodynamics", config_use_column_vertical_thermodynamics)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nVerticesSolve", nVerticesSolve)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
    call MPAS_pool_get_array(mesh, "fVertex", fVertex)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

    call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)
    call MPAS_pool_get_array(tracers_aggregate, "iceVolumeCell", iceVolumeCell)
    call MPAS_pool_get_array(tracers_aggregate, "snowVolumeCell", snowVolumeCell)

    call MPAS_pool_get_array(icestate, "iceAreaVertex", iceAreaVertex)   
    call MPAS_pool_get_array(icestate, "totalMassCell", totalMassCell)
    call MPAS_pool_get_array(icestate, "iceAreaCellInitial", iceAreaCellInitial)
    call MPAS_pool_get_array(icestate, "totalMassVertex", totalMassVertex)
    call MPAS_pool_get_array(icestate, "openWaterArea", openWaterArea)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)
    call MPAS_pool_get_array(velocity_solver, "uVelocityInitial", uVelocityInitial)
    call MPAS_pool_get_array(velocity_solver, "vVelocityInitial", vVelocityInitial)
    call MPAS_pool_get_array(velocity_solver, "solveStress", solveStress)
    call MPAS_pool_get_array(velocity_solver, "solveVelocity", solveVelocity)
    call MPAS_pool_get_array(velocity_solver, "solveVelocityPrevious", solveVelocityPrevious)    
    call MPAS_pool_get_array(velocity_solver, "icePressure", icePressure)
    call MPAS_pool_get_array(velocity_solver, "airStressCellU", airStressCellU)
    call MPAS_pool_get_array(velocity_solver, "airStressCellV", airStressCellV)
    call MPAS_pool_get_array(velocity_solver, "airStressVertexU", airStressVertexU)
    call MPAS_pool_get_array(velocity_solver, "airStressVertexV", airStressVertexV)
    call MPAS_pool_get_array(velocity_solver, "uOceanVelocityVertex", uOceanVelocityVertex)
    call MPAS_pool_get_array(velocity_solver, "vOceanVelocityVertex", vOceanVelocityVertex)
    call MPAS_pool_get_array(velocity_solver, "totalMassVertexfVertex", totalMassVertexfVertex)
    call MPAS_pool_get_array(velocity_solver, "oceanStressU", oceanStressU)
    call MPAS_pool_get_array(velocity_solver, "oceanStressV", oceanStressV)
    call MPAS_pool_get_array(velocity_solver, "surfaceTiltForceU", surfaceTiltForceU)
    call MPAS_pool_get_array(velocity_solver, "surfaceTiltForceV", surfaceTiltForceV)
    call MPAS_pool_get_array(velocity_solver, "stressDivergenceU", stressDivergenceU)
    call MPAS_pool_get_array(velocity_solver, "stressDivergenceV", stressDivergenceV)

    call MPAS_pool_get_array(atmos_coupling, "uAirVelocity", uAirVelocity)
    call MPAS_pool_get_array(atmos_coupling, "vAirVelocity", vAirVelocity)
    call MPAS_pool_get_array(atmos_coupling, "airDensity", airDensity)

    call MPAS_pool_get_array(ocean_coupling, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(ocean_coupling, "vOceanVelocity", vOceanVelocity)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    ! aggregate categories for area and volume
    call aggregate_mass_and_area(&
         nCells, &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         totalMassCell)

    ! interpolate area and mass from cells to vertices
    call cice_interpolate_cell_to_vertex(&
         mesh, &
         iceAreaVertex, &
         iceAreaCellInitial)

    call cice_interpolate_cell_to_vertex(&
         mesh, &
         totalMassVertex, &
         totalMassCell)

    ! calculate computational masks
    call stress_calculation_mask(&
         nCells, &
         mesh, &
         solveStress, &
         iceAreaCellInitial, &
         totalMassCell)

    call velocity_calculation_mask(&
         nVerticesSolve, &
         mesh, &
         solveVelocity,   &
         iceAreaVertex,   &
         totalMassVertex, &
         interiorVertex)

    ! interpolate ocean velocities
    call cice_interpolate_cell_to_vertex(&
         mesh, &
         uOceanVelocityVertex, &
         uOceanVelocity)

    call cice_interpolate_cell_to_vertex(&
         mesh, &
         vOceanVelocityVertex, &
         vOceanVelocity)
    
    ! set new ice velocities
    call new_ice_velocities(&
         nVerticesSolve, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         solveVelocity, &
         solveVelocityPrevious, &
         uVelocity, &
         vVelocity, &
         uVelocityInitial, &
         vVelocityInitial, &
         stressDivergenceU, &
         stressDivergenceV, &
         oceanStressU, &
         oceanStressV)

    ! calculate the ice strength
    call cice_ice_strength(&
         nCells, &     
         nCategories, &
         icePressure, &
         iceAreaCell, &
         iceVolumeCell, &
         openWaterArea, &
         iceAreaCategory, &
         iceVolumeCategory, &
         solveStress, &
         config_use_column_package, &
         config_use_column_vertical_thermodynamics)

    ! calculate the air stress
    if (.not. config_use_column_package .or. &
       (config_use_column_package .and. .not. config_use_column_vertical_thermodynamics)) then
       call cice_air_stress(&
            nCells, &
            airStressCellU, &
            airStressCellV, &
            uAirVelocity, &
            vAirVelocity, &
            airDensity, &
            iceAreaCell)
    endif

    !airStressCellU = 0.0_RKIND
    !airStressCellV = 0.0_RKIND

    call cice_interpolate_cell_to_vertex(&
         mesh, &
         airStressVertexU, &
         airStressCellU)

    call cice_interpolate_cell_to_vertex(&
         mesh, &
         airStressVertexV, &
         airStressCellV)

    ! calculate the coriolis force coefficient
    call cice_coriolis_force_coefficient(&
         nVerticesSolve, &
         totalMassVertexfVertex, &
         totalMassVertex, &
         fVertex)

    ! calculate the ocean stress
    call cice_ocean_stress(& 
         nVerticesSolve, &
         oceanStressU, &
         oceanStressV, &
         uOceanVelocityVertex, & 
         vOceanVelocityVertex, &
         fVertex, &
         solveVelocity)

    ! calculate the surface tilt force
    call cice_surface_tilt(&
         nVerticesSolve, &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         uOceanVelocityVertex, & 
         vOceanVelocityVertex, &
         totalMassVertex, &
         fVertex, &
         solveVelocity)

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

  subroutine velocity_solver_post_subcycle(block)

    use cice_velocity_solver_shared, only: &
         cice_ocean_stress_final

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         boundary, &
         configs, &
         icestate, &
         velocity_weak, &
         velocity_variational, &
         velocity_solver

    integer, pointer :: &
         nVerticesSolve

    integer, dimension(:), pointer :: &
         solveVelocity

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaVertex, &
         replacementPressure, &
         principalStress1, &
         principalStress2, &
         uVelocity, &
         vVelocity, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         oceanStressU, &
         oceanStressV, &
         oceanStressCoeff, &
         oceanStressCellU, &
         oceanStressCellV, &
         fVertex

    configs => block % configs

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
    call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocity_weak)
    call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocity_variational)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)

    call MPAS_pool_get_dimension(mesh, "nVerticesSolve", nVerticesSolve)

    call MPAS_pool_get_array(mesh, "fVertex", fVertex)

    call MPAS_pool_get_array(icestate, "iceAreaVertex", iceAreaVertex)

    call MPAS_pool_get_array(velocity_solver, "replacementPressure", replacementPressure)
    call MPAS_pool_get_array(velocity_solver, "principalStress1", principalStress1)
    call MPAS_pool_get_array(velocity_solver, "principalStress2", principalStress2)
    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)
    call MPAS_pool_get_array(velocity_solver, "uOceanVelocityVertex", uOceanVelocityVertex)
    call MPAS_pool_get_array(velocity_solver, "vOceanVelocityVertex", vOceanVelocityVertex)
    call MPAS_pool_get_array(velocity_solver, "oceanStressU", oceanStressU)
    call MPAS_pool_get_array(velocity_solver, "oceanStressV", oceanStressV)
    call MPAS_pool_get_array(velocity_solver, "oceanStressCoeff", oceanStressCoeff)
    call MPAS_pool_get_array(velocity_solver, "oceanStressCellU", oceanStressCellU)
    call MPAS_pool_get_array(velocity_solver, "oceanStressCellV", oceanStressCellV)
    call MPAS_pool_get_array(velocity_solver, "solveVelocity", solveVelocity)

    call final_divergence_shear(&
         block)
    
    call principal_stresses_driver(&
         mesh, &
         velocity_weak, &
         velocity_variational, &
         replacementPressure, &
         configs, &
         principalStress1, &
         principalStress2)

    call cice_ocean_stress_final(&
         mesh, &
         boundary, &
         nVerticesSolve, &
         iceAreaVertex, &
         fVertex, &
         uVelocity, &
         vVelocity, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         oceanStressU, &
         oceanStressV, &
         oceanStressCoeff, &
         oceanStressCellU, &
         oceanStressCellV, &
         solveVelocity)

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
  
  subroutine final_divergence_shear(&
       block)

    use cice_velocity_solver_variational, only: &
         cice_final_divergence_shear_variational

    use cice_velocity_solver_weak, only: &
         cice_final_divergence_shear_weak

    type(block_type), intent(inout) :: &
         block
    
    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    call MPAS_pool_get_config(block % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)
    
    if (trim(config_stress_divergence_scheme) == "weak") then

       call cice_final_divergence_shear_weak(block)
       
    else if (trim(config_stress_divergence_scheme) == "variational") then

       call cice_final_divergence_shear_variational(block)

    endif

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

  subroutine principal_stresses_driver(&
       mesh, &
       velocity_weak, &
       velocity_variational, &
       replacementPressure, &
       configs, &
       principalStress1, &
       principalStress2)!{{{

    use cice_velocity_solver_shared, only: &
         cice_principal_stresses         

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh, &                 !< Input: 
         velocity_weak, &        !< Input: 
         velocity_variational, & !< Input: 
         configs                 !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         replacementPressure !< Input: 
    
    real(kind=RKIND), dimension(:), intent(out) :: &
         principalStress1, & !< Output: 
         principalStress2    !< Output: 

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

    call MPAS_pool_get_config(configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)

    ! calculate the principal stresses
    if (trim(config_stress_divergence_scheme) == "weak") then

       call MPAS_pool_get_array(velocity_weak, "stress11", stress11Weak) 
       call MPAS_pool_get_array(velocity_weak, "stress22", stress22Weak) 
       call MPAS_pool_get_array(velocity_weak, "stress12", stress12Weak) 

       call cice_principal_stresses(&
            nCellsSolve, &
            principalStress1, &
            principalStress2, &
            stress11Weak, &
            stress22Weak, &
            stress12Weak, &
            replacementPressure)

    else if (trim(config_stress_divergence_scheme) == "variational") then

       call MPAS_pool_get_array(velocity_variational, "stress11", stress11Var) 
       call MPAS_pool_get_array(velocity_variational, "stress22", stress22Var)
       call MPAS_pool_get_array(velocity_variational, "stress12", stress12Var)

       call cice_principal_stresses(&
            nCellsSolve, &
            principalStress1, &
            principalStress2, &
            stress11Var(1,:), &
            stress22Var(1,:), &
            stress12Var(1,:), &
            replacementPressure)

    endif

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

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         icestate, &
         boundary, &
         velocity_solver, &
         velocity_weak, &
         velocity_variational

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme, &
         config_variational_basis

    logical, pointer :: &
         config_revised_evp

    integer, pointer :: &
         nVerticesSolve

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveStress

    type(field1DReal), pointer :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocityInitial, &
         vVelocityInitial, &
         iceAreaVertex, &
         totalMassVertex, &
         fVertex, &
         icePressure, &
         replacementPressure, &
         stressDivergenceU, &
         stressDivergenceV, &
         airStressVertexU, &
         airStressVertexV, &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         totalMassVertexfVertex, &
         oceanStressU, &
         oceanStressV, &
         oceanStressCoeff, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex
         
    real(kind=RKIND), pointer :: &
         elasticTimeStep, &
         config_dt

    type(MPAS_pool_type), pointer :: &
         velocity_solver_fields

    type(field1DReal), pointer :: &
         uVelocityField, &
         vvelocityField

    call MPAS_pool_get_config(domain % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)
    call MPAS_pool_get_config(domain % configs, "config_variational_basis", config_variational_basis)
    call MPAS_pool_get_config(domain % configs, "config_revised_evp", config_revised_evp)
    call MPAS_pool_get_config(domain % configs, "config_dt", config_dt)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
       call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
       call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocity_weak)
       call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocity_variational)

       call MPAS_pool_get_array(icestate, "iceAreaVertex", iceAreaVertex)
       call MPAS_pool_get_array(icestate, "totalMassVertex", totalMassVertex)

       call MPAS_pool_get_dimension(mesh, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_array(velocity_solver, "elasticTimeStep", elasticTimeStep)
       call MPAS_pool_get_array(velocity_solver, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocity_solver, "solveStress", solveStress)
       call MPAS_pool_get_field(velocity_solver, "uVelocity", uVelocity)
       call MPAS_pool_get_field(velocity_solver, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocity_solver, "uVelocityInitial", uVelocityInitial)
       call MPAS_pool_get_array(velocity_solver, "vVelocityInitial", vVelocityInitial)
       call MPAS_pool_get_array(velocity_solver, "icePressure", icePressure)
       call MPAS_pool_get_array(velocity_solver, "replacementPressure", replacementPressure)
       call MPAS_pool_get_array(velocity_solver, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocity_solver, "stressDivergenceV", stressDivergenceV)
       call MPAS_pool_get_array(velocity_solver, "airStressVertexU", airStressVertexU)
       call MPAS_pool_get_array(velocity_solver, "airStressVertexV", airStressVertexV)
       call MPAS_pool_get_array(velocity_solver, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocity_solver, "surfaceTiltForceV", surfaceTiltForceV)
       call MPAS_pool_get_array(velocity_solver, "totalMassVertexfVertex", totalMassVertexfVertex)
       call MPAS_pool_get_array(velocity_solver, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocity_solver, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocity_solver, "oceanStressCoeff", oceanStressCoeff)
       call MPAS_pool_get_array(velocity_solver, "uOceanVelocityVertex", uOceanVelocityVertex)
       call MPAS_pool_get_array(velocity_solver, "vOceanVelocityVertex", vOceanVelocityVertex)
       
       if (trim(config_stress_divergence_scheme) == "weak") then

          call cice_internal_stress_weak(&
               mesh, &
               velocity_weak, &
               uVelocity % array, &
               vVelocity % array, &
               icePressure, &
               replacementPressure, &
               elasticTimeStep, &
               config_dt, &
               config_revised_evp, &
               stressDivergenceU, &
               stressDivergenceV, &
               solveStress, &
               solveVelocity)

       else if (trim(config_stress_divergence_scheme) == "variational") then

          call cice_internal_stress_variational(&
               mesh, &
               velocity_variational, &
               uVelocity % array, &
               vVelocity % array, &
               icePressure, &
               replacementPressure, &
               elasticTimeStep, &
               config_dt, &
               config_revised_evp, &
               stressDivergenceU, &
               stressDivergenceV, &
               solveStress, &
               solveVelocity)

       endif

       !stressDivergenceU = 0.0_RKIND
       !stressDivergenceV = 0.0_RKIND

       call cice_ocean_stress_coefficient(& 
            nVerticesSolve, &
            oceanStressCoeff, &
            uOceanVelocityVertex, & 
            vOceanVelocityVertex, &
            uVelocity % array, &   
            vVelocity % array, &
            iceAreaVertex, &
            solveVelocity)

       if (.not. config_revised_evp) then
       
          call cice_solve_velocity(&
               nVerticesSolve, &
               solveVelocity, &
               uVelocity % array, &   
               vVelocity % array, &
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
               oceanStressCoeff, &
               elasticTimeStep, &
               mesh)

       else

          call cice_solve_velocity_revised(&
               nVerticesSolve, &
               solveVelocity, &
               uVelocity % array, &   
               vVelocity % array, &
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
               oceanStressCoeff, &
               config_dt, &
               mesh)

       endif

       block => block % next
    end do
      
    ! halo exchange
    call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocity_solver_fields)

    call MPAS_pool_get_field(velocity_solver_fields, "uVelocity", uVelocityField)
    call MPAS_pool_get_field(velocity_solver_fields, "vVelocity", vVelocityField)

    call MPAS_dmpar_exch_halo_field(uVelocityField)
    call MPAS_dmpar_exch_halo_field(vVelocityField)

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

  subroutine aggregate_mass_and_area(&
       nPoints, &
       iceAreaCategory, &
       iceVolumeCategory, &
       snowVolumeCategory, &
       iceAreaCell, &
       iceVolumeCell, &
       snowVolumeCell, &
       totalMassCell)!{{{

    use cice_constants, only: &
         ciceDensityIce, &
         ciceDensitySnow

    integer, intent(in) :: &
         nPoints

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         iceAreaCategory, &   !< Input: 
         iceVolumeCategory, & !< Input: 
         snowVolumeCategory   !< Input: 

    real(kind=RKIND), dimension(:), intent(out) :: &
         totalMassCell, & !< Output: 
         iceAreaCell, &   !< Output: 
         iceVolumeCell, & !< Output: 
         snowVolumeCell   !< Output: 

    integer :: &
         iPoint

    do iPoint = 1, nPoints

       iceAreaCell(iPoint)    = sum(iceAreaCategory(1,:,iPoint))
       iceVolumeCell(iPoint)  = sum(iceVolumeCategory(1,:,iPoint))
       snowVolumeCell(iPoint) = sum(snowVolumeCategory(1,:,iPoint))

       totalMassCell(iPoint)  = iceVolumeCell(iPoint)  * ciceDensityIce + &
                                snowVolumeCell(iPoint) * ciceDensitySnow

    enddo ! nPoints

  end subroutine aggregate_mass_and_area!}}}

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

  subroutine stress_calculation_mask(&
       nPoints, &
       mesh, &
       solveStress, &
       areaCell, &
       massCell)!{{{

    use cice_velocity_solver_shared, only: &
         ciceAreaMinimum, &
         ciceMassMinimum

    integer, intent(in) :: &
         nPoints !< Input: 

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    integer, dimension(:), intent(out) :: &
         solveStress !< Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         areaCell, & !< Input: 
         massCell    !< Input: 

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

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "cellsOnCell", cellsOnCell)

    do iCell = 1, nPoints

       solveStress(iCell) = 0

       if (areaCell(iCell) > ciceAreaMinimum .and. &
           massCell(iCell) > ciceMassMinimum) then

          ! this cell has sufficient ice
          solveStress(iCell) = 1

       else

          ! test neighbouring cells to see if have sufficient ice
          do iCellOnCell = 1, nEdgesOnCell(iCell)

             iCellNeighbour = cellsOnCell(iCellOnCell,iCell)

             if (areaCell(iCellNeighbour) > ciceAreaMinimum .and. &
                 massCell(iCellNeighbour) > ciceMassMinimum) then

                solveStress(iCell) = 1
                exit

             endif

          enddo ! iCellOnCell

       endif

    enddo ! iCell

    do iCell = nPoints+1, nCells

       solveStress(iCell) = 0

    enddo ! iCell

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

  subroutine velocity_calculation_mask(&
       nPoints, &
       mesh, &
       solveVelocity, &
       areaVertex, &
       massVertex, &
       interiorVertex)!{{{

    use cice_velocity_solver_shared, only: &
         ciceAreaMinimum, &
         ciceMassMinimum

    integer, intent(in) :: &
         nPoints

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    integer, dimension(:), intent(out) :: &
         solveVelocity !< Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         areaVertex, & !< Input: 
         massVertex    !< Input: 

    integer, dimension(:), intent(in) :: &
         interiorVertex !< Input: 

    integer :: &
         iVertex

    integer, pointer :: &
         nVertices

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    do iVertex = 1, nPoints

       solveVelocity(iVertex) = 0

       if (interiorVertex(iVertex) == 1 .and. &
           areaVertex(iVertex) > ciceAreaMinimum .and. &
           massVertex(iVertex) > ciceMassMinimum) then

          ! this vertex has sufficient ice
          solveVelocity(iVertex) = 1

       endif

    enddo ! iVertex

    do iVertex = nPoints+1, nVertices

       solveVelocity(iVertex) = 0

    enddo ! iVertex

  end subroutine velocity_calculation_mask!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 29th June 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine new_ice_velocities(&
       nVertices, &
       uOceanVelocityVertex, &
       vOceanVelocityVertex, &
       solveVelocity, &
       solveVelocityPrevious, &
       uVelocity, &
       vVelocity, &
       uVelocityInitial, &
       vVelocityInitial, &
       stressDivergenceU, &
       stressDivergenceV, &
       oceanStressU, &
       oceanStressV)!{{{

    integer, intent(in) :: &
         nVertices
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         uOceanVelocityVertex, &
         vOceanVelocityVertex

    integer, dimension(:), intent(in) :: &
         solveVelocity

    integer, dimension(:), intent(inout) :: &
         solveVelocityPrevious

    real(kind=RKIND), dimension(:), intent(inout) :: &
         uVelocity, &
         vVelocity, &
         uVelocityInitial, &
         vVelocityInitial, &
         stressDivergenceU, &
         stressDivergenceV, &
         oceanStressU, &
         oceanStressV
         
    integer :: &
         iVertex

    do iVertex = 1, nVertices

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
