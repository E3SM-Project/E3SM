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

  use mpas_grid_types
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
       block, &
       dt)!{{{

    use cice_velocity_solver_shared, only: &
         cice_init_evp

    use cice_velocity_solver_weak, only: &
         cice_init_velocity_solver_weak

    use cice_velocity_solver_variational, only: &
         cice_init_velocity_solver_variational

    type (block_type), intent(inout) :: &
         block !< Input/Output: 

    real(kind=RKIND), intent(in) :: &
         dt !< Input: 

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme, &
         config_variational_basis

    logical, pointer :: &
         config_rotate_cartesian_grid, &
         config_include_metric_terms

    type (MPAS_pool_type), pointer :: &
         mesh, &
         boundary, &
         velocity_weak, &
         velocity_variational, &
         velocity_pwl, &
         configs

    configs => block % configs

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocity_weak)
    call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocity_variational)

    call MPAS_pool_get_config(configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)
    call MPAS_pool_get_config(configs, "config_variational_basis", config_variational_basis)
    call MPAS_pool_get_config(configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)
    call MPAS_pool_get_config(configs, "config_include_metric_terms", config_include_metric_terms)
    
    call cice_init_evp(dt, configs)

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
       block, &
       clock, &
       dt)!{{{

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    use cice_velocity_solver_shared, only: &
         cice_ice_strength, &           
         cice_air_stress, &
         cice_ocean_stress, &
         cice_ocean_stress_final, &
         cice_surface_tilt, &
         cice_coriolis_force_coefficient, &
         ciceDynamicsTimeStep, &
         ciceElasticTimeStep

    use mpas_constants, only: &
         pii

    use cice_debug, only: &
         cice_time_output_variable_real, &
         cice_time_output_variable_int

    type (block_type), intent(inout) :: &
         block !< Input/Output: 

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input: 

    real(kind=RKIND), intent(in) :: &
         dt !< Input: 

    type (MPAS_pool_type), pointer :: &
         mesh, &
         configs, &
         tracers, &
         icestate, &
         boundary, &
         atmos_forcing, &
         ocean_forcing, &
         velocity_solver, &
         velocity_weak, &
         velocity_variational

    integer, pointer :: &
         nCells, &
         nCellsSolve, &
         nVertices, &
         nVerticesSolve

    integer, dimension(:), pointer :: &
         solveStress, &
         solveVelocity, &
         interiorVertex

    type(field3DReal), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         totalMassCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         iceAreaVertex, &
         totalMassVertex, &
         icePressure, &
         airStressVertexU, &
         airStressVertexV, &
         oceanStressU, &
         oceanStressV, &
         oceanStressCoeff, &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         totalMassVertexfVertex, &
         principalStress1, &
         principalStress2, &
         replacementPressure, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         airDensity, &
         uAirVelocity, &
         vAirVelocity, &
         uVelocity, &
         vVelocity, &
         fVertex, &
         latCell, &
         airStressCellU, &
         airStressCellV, &
         uOceanVelocity, &
         vOceanVelocity

    configs => block % configs

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmos_forcing)
    call MPAS_pool_get_subpool(block % structs, "ocean_forcing", ocean_forcing)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
    call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocity_weak)
    call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocity_variational)

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nVerticesSolve", nVerticesSolve)
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_field(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_field(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_field(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

    call MPAS_pool_get_array(icestate, "iceAreaCell", iceAreaCell)
    call MPAS_pool_get_array(icestate, "totalMassCell", totalMassCell)

    call MPAS_pool_get_array(icestate, "iceVolumeCell", iceVolumeCell)
    call MPAS_pool_get_array(icestate, "snowVolumeCell", snowVolumeCell)
    call MPAS_pool_get_array(icestate, "iceAreaVertex", iceAreaVertex)
    call MPAS_pool_get_array(icestate, "totalMassVertex", totalMassVertex)

    call MPAS_pool_get_array(velocity_solver, "solveStress", solveStress)
    call MPAS_pool_get_array(velocity_solver, "solveVelocity", solveVelocity)
    call MPAS_pool_get_array(velocity_solver, "icePressure", icePressure)
    call MPAS_pool_get_array(velocity_solver, "airStressCellU", airStressCellU)
    call MPAS_pool_get_array(velocity_solver, "airStressCellV", airStressCellV)
    call MPAS_pool_get_array(velocity_solver, "airStressVertexU", airStressVertexU)
    call MPAS_pool_get_array(velocity_solver, "airStressVertexV", airStressVertexV)
    call MPAS_pool_get_array(velocity_solver, "oceanStressU", oceanStressU)
    call MPAS_pool_get_array(velocity_solver, "oceanStressV", oceanStressV)
    call MPAS_pool_get_array(velocity_solver, "oceanStressCoeff", oceanStressCoeff)
    call MPAS_pool_get_array(velocity_solver, "surfaceTiltForceU", surfaceTiltForceU)
    call MPAS_pool_get_array(velocity_solver, "surfaceTiltForceV", surfaceTiltForceV)
    call MPAS_pool_get_array(velocity_solver, "totalMassVertexfVertex", totalMassVertexfVertex)
    call MPAS_pool_get_array(velocity_solver, "principalStress1", principalStress1)
    call MPAS_pool_get_array(velocity_solver, "principalStress2", principalStress2)
    call MPAS_pool_get_array(velocity_solver, "replacementPressure", replacementPressure)
    call MPAS_pool_get_array(velocity_solver, "uOceanVelocityVertex", uOceanVelocityVertex)
    call MPAS_pool_get_array(velocity_solver, "vOceanVelocityVertex", vOceanVelocityVertex)
    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)

    call MPAS_pool_get_array(ocean_forcing, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(ocean_forcing, "vOceanVelocity", vOceanVelocity)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call MPAS_pool_get_array(atmos_forcing, "uAirVelocity", uAirVelocity)
    call MPAS_pool_get_array(atmos_forcing, "vAirVelocity", vAirVelocity)
    call MPAS_pool_get_array(atmos_forcing, "airDensity", airDensity)

    call MPAS_pool_get_array(mesh, "fVertex", fVertex)

    !call cice_time_output_variable_real(mesh, uVelocity, clock, "uVelocity0")

    ! halo updates
    call mpas_dmpar_exch_halo_field(iceAreaCategory)
    call mpas_dmpar_exch_halo_field(iceVolumeCategory)
    call mpas_dmpar_exch_halo_field(snowVolumeCategory)

    ! aggregate categories for area and volume
    call aggregate_mass_and_area(&
         nCells, &
         iceAreaCategory % array, &
         iceVolumeCategory % array, &
         snowVolumeCategory % array, &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         totalMassCell)

    ! interpolate area and mass from cells to vertices
    call cice_interpolate_cell_to_vertex(&
         mesh, &
         iceAreaVertex, &
         iceAreaCell)

    call cice_interpolate_cell_to_vertex(&
         mesh, &
         totalMassVertex, &
         totalMassCell)

    ! calculate computational masks
    call stress_calculation_mask(&
         nCells, &
         mesh, &
         solveStress, &
         iceAreaCell, &
         totalMassCell)

    call velocity_calculation_mask(&
         nVerticesSolve, &
         mesh, &
         solveVelocity,   &
         iceAreaVertex,   &
         totalMassVertex, &
         interiorVertex)

    ! calculate the ice strength
    call cice_ice_strength(&
         nCells, &     
         icePressure, &
         iceAreaCell, &
         iceVolumeCell, &
         solveStress)

    ! calculate the air stress
    call cice_air_stress(&
         nCells, &
         airStressCellU, &
         airStressCellV, &
         uAirVelocity, & 
         vAirVelocity, &
         airDensity, &
         iceAreaCell)

    call cice_interpolate_cell_to_vertex(&
         mesh, &
         airStressVertexU, &
         airStressCellU)

    call cice_interpolate_cell_to_vertex(&
         mesh, &
         airStressVertexV, &
         airStressCellV)

    call cice_interpolate_cell_to_vertex(&
         mesh, &
         uOceanVelocityVertex, &
         uOceanVelocity)

    call cice_interpolate_cell_to_vertex(&
         mesh, &
         vOceanVelocityVertex, &
         vOceanVelocity)

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

    ! subcycle the dynamics
    call subcycle_velocity_solver(&
         mesh, &
         clock, &
         icestate, &
         boundary, &
         velocity_solver, &
         velocity_weak, &
         velocity_variational, &
         ciceDynamicsTimeStep, &
         ciceElasticTimeStep, &
         configs)

    call principal_stresses_driver(&
         mesh, &
         velocity_weak, &
         velocity_variational, &
         replacementPressure, &
         configs, &
         principalStress1, &
         principalStress2)

    call cice_ocean_stress_final(&
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
         solveVelocity)

    !call cice_time_output_variable_real(mesh, uVelocity, clock, "uVelocity")
    !call cice_time_output_variable_real(mesh, totalMassVertex, clock, "totalMassVertex")
    !call cice_time_output_variable_real(mesh, iceAreaVertex, clock, "iceAreaVertex")
    !call cice_time_output_variable_real(mesh, oceanStressU, clock, "oceanStressU")
    !call cice_time_output_variable_real(mesh, airStressCellU, clock, "airStressCellU")
    !call cice_time_output_variable_real(mesh, icePressure, clock, "icePressure")

    !call cice_time_output_variable_int(mesh, solveVelocity, clock, "solveVelocity")
    !call cice_time_output_variable_int(mesh, solveStress, clock, "solveStress")

  end subroutine cice_run_velocity_solver!}}}

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
       mesh, &
       clock, &
       icestate, &
       boundary,  &
       velocity_solver, &
       velocity_weak, &
       velocity_variational, &
       dtDynamics, &
       dtElastic, &
       configs)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh, & !< Input: 
         configs !< Input: 

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input: 

    type(MPAS_pool_type), pointer :: &
         icestate, &          !< Input/Output: 
         boundary, &          !< Input/Output: 
         velocity_solver, &   !< Input/Output: 
         velocity_weak, &     !< Input/Output: 
         velocity_variational !< Input/Output: 

    real(kind=RKIND), intent(in) :: &
         dtDynamics, & !< Input: 
         dtElastic     !< Input: 
 
    integer, pointer :: &
         config_elastic_subcycle_number

    integer :: &
         iElasticSubcycle

    call MPAS_pool_get_config(configs, "config_elastic_subcycle_number", config_elastic_subcycle_number)

    do iElasticSubcycle = 1, config_elastic_subcycle_number

       call single_subcycle_velocity_solver(&
            mesh, &
            clock, &
            icestate, &
            boundary, &
            velocity_solver, &
            velocity_weak, &
            velocity_variational, &
            dtDynamics, &
            dtElastic, &
            configs, &
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
       mesh, &
       clock, &
       icestate, &
       boundary,  &
       velocity_solver, &
       velocity_weak, &
       velocity_variational, &
       dtDynamics, &
       dtElastic, &
       configs, &
       iElasticSubcycle)!{{{                               

    use cice_velocity_solver_shared, only: &
         cice_ocean_stress_coefficient, &
         cice_solve_velocity

    use cice_velocity_solver_weak, only: &
         cice_internal_stress_weak

    use cice_velocity_solver_variational, only: &
         cice_internal_stress_variational

    use cice_debug, only: &
         cice_time_output_variable_real

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh, & !< Input: 
         configs !< Input: 

    type(MPAS_Clock_type), intent(in) :: &
         clock

    type(MPAS_pool_type), pointer :: &
         icestate, &          !< Input/Output: 
         boundary, &          !< Input/Output: 
         velocity_solver, &   !< Input/Output: 
         velocity_weak, &     !< Input/Output: 
         velocity_variational !< Input/Output: 

    real(kind=RKIND), intent(in) :: &
         dtDynamics, & !< Input: 
         dtElastic     !< Input: 

    integer, intent(in) :: &
         iElasticSubcycle !< Input: !! testing

    integer :: &
         iVertex !! testing

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme, &
         config_variational_basis

    logical, pointer :: &
         config_evp_damping

    integer, pointer :: &
         nVerticesSolve

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveStress

    type(field1DReal), pointer :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:), pointer :: &
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

    call MPAS_pool_get_config(configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)
    call MPAS_pool_get_config(configs, "config_variational_basis", config_variational_basis)
    call MPAS_pool_get_config(configs, "config_evp_damping", config_evp_damping)

    call MPAS_pool_get_array(icestate, "iceAreaVertex", iceAreaVertex)
    call MPAS_pool_get_array(icestate, "totalMassVertex", totalMassVertex)

    call MPAS_pool_get_dimension(mesh, "nVerticesSolve", nVerticesSolve)

    call MPAS_pool_get_array(velocity_solver, "solveVelocity", solveVelocity)
    call MPAS_pool_get_array(velocity_solver, "solveStress", solveStress)
    call MPAS_pool_get_field(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_field(velocity_solver, "vVelocity", vVelocity)
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
            dtElastic, &
            config_evp_damping, &
            stressDivergenceU, &
            stressDivergenceV, &
            solveStress, &
            solveVelocity)

    else if (trim(config_variational_basis) == "wachspress") then

       call cice_internal_stress_variational(&
            mesh, &
            velocity_variational, &
            uVelocity % array, &
            vVelocity % array, &
            icePressure, &
            replacementPressure, &
            dtElastic, &
            config_evp_damping, &
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

    call cice_solve_velocity(&
         nVerticesSolve, &
         solveVelocity, &
         uVelocity % array, &   
         vVelocity % array, &
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
         dtElastic, &
         mesh)

    !call cice_time_output_variable_real(mesh, stressDivergenceU, clock, "stressDivergenceU")

    ! halo exchange
    call mpas_dmpar_exch_halo_field(uVelocity)
    call mpas_dmpar_exch_halo_field(vVelocity)

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

!-----------------------------------------------------------------------

end module cice_velocity_solver
