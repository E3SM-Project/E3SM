module cice_dynamics_hex

  use mpas_grid_types

  implicit none

  private
  public :: &
       init_dynamics_hex, &
       run_dynamics_hex
  
contains

  !-------------------------------------------------------------
  ! Initialization
  !-------------------------------------------------------------
  
  subroutine init_dynamics_hex(mesh,     dt,     &
                               boundary,         &
                               hexweak,  hexvar, &
                               hexwach,  hexpwl, &
                               configs)

    use cice_dynamics_shared, only: init_evp
    use cice_dyn_hex_weak, only: init_dynamics_hex_weak
    use cice_dyn_hex_wachspress, only: init_dynamics_hex_wachspress
    use cice_dyn_hex_pwl, only: init_dynamics_hex_pwl

    use cice_diagnostics, only: plot_pop_grid_vertices

    type(MPAS_pool_type), pointer, intent(inout) :: mesh
    type(MPAS_pool_type), pointer, intent(in) :: configs

    type(MPAS_pool_type), pointer :: hexweak
    type(MPAS_pool_type), pointer :: hexvar
    type(MPAS_pool_type), pointer :: hexwach
    type(MPAS_pool_type), pointer :: hexpwl
    type(MPAS_pool_type), pointer :: boundary

    real(kind=RKIND), intent(in) :: dt

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme, &
         config_variational_basis

    logical, pointer :: &
         config_rotate_cartesian_grid

    call MPAS_pool_get_config(configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)
    call MPAS_pool_get_config(configs, "config_variational_basis", config_variational_basis)
    call MPAS_pool_get_config(configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)
    
    call init_evp(dt, configs)

    if (trim(config_stress_divergence_scheme) == "weak") then

       call init_dynamics_hex_weak(mesh, boundary, hexweak, config_rotate_cartesian_grid)

    else if (trim(config_stress_divergence_scheme) == "variational" .and. &
             trim(config_variational_basis) == "wachspress") then

       call init_dynamics_hex_wachspress(mesh, hexvar, hexwach, boundary)

    else if (trim(config_stress_divergence_scheme) == "variational" .and. &
             trim(config_variational_basis) == "pwl") then

       call init_dynamics_hex_pwl(mesh, hexvar, hexpwl, boundary)

    endif

  end subroutine init_dynamics_hex

  !-------------------------------------------------------------
  ! Time stepping
  !-------------------------------------------------------------

  subroutine run_dynamics_hex(mesh, &
                              tracers,       icestate,      &
                              boundary,                     &
                              atmos_forcing, ocean_forcing, &
                              hexdyn,                       &
                              hexweak,       hexvar,        &
                              hexwach,       hexpwl,        &
                              dt,            configs,       &
                              nstep)

    use cice_mesh, only: &
         interpolate_cell_to_vertex

    use cice_dynamics_shared, only: &
         ice_strength, &           
         air_stress, &
         ocean_stress, &
         ocean_stress_final, &
         surface_tilt, &
         coriolis_force_coefficient, &
         dtDynamics, &
         dtElastic, &
         elapsedTime
    
    use cice_diagnostics, only: &
         plot_vertices_planar         

    use mpas_constants, only: pii

    type(MPAS_pool_type), pointer, intent(in) :: mesh
    type(MPAS_pool_type), pointer, intent(in) :: configs
    type(MPAS_pool_type), pointer :: tracers
    type(MPAS_pool_type), pointer :: icestate
    type(MPAS_pool_type), pointer :: boundary
    type(MPAS_pool_type), pointer :: atmos_forcing 
    type(MPAS_pool_type), pointer :: ocean_forcing    
    type(MPAS_pool_type), pointer :: hexdyn
    type(MPAS_pool_type), pointer :: hexweak
    type(MPAS_pool_type), pointer :: hexvar
    type(MPAS_pool_type), pointer :: hexwach
    type(MPAS_pool_type), pointer :: hexpwl

    real(kind=RKIND), intent(in) :: &
         dt

    integer, intent(in) :: &
         nstep

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         solveStress, &
         solveVelocity, &
         interiorVertex

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         iceAreaVertex, &
         totalMassVertex, &
         totalMassCell, &
         icePressure, &
         airStressCellU, &
         airStressCellV, &
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
         uOceanVelocity, &
         vOceanVelocity, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         airDensity, &
         uAirVelocity, &
         vAirVelocity, &
         uVelocity, &
         vVelocity, &
         fVertex, &
         latCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    call MPAS_pool_set_error_level(MPAS_POOL_FATAL) ! SILENT, WARN and FATAL

    ! init variables
    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

    call MPAS_pool_get_array(icestate, "iceAreaCell", iceAreaCell)
    call MPAS_pool_get_array(icestate, "iceVolumeCell", iceVolumeCell)
    call MPAS_pool_get_array(icestate, "snowVolumeCell", snowVolumeCell)
    call MPAS_pool_get_array(icestate, "iceAreaVertex", iceAreaVertex)
    call MPAS_pool_get_array(icestate, "totalMassVertex", totalMassVertex)
    call MPAS_pool_get_array(icestate, "totalMassCell", totalMassCell)

    call MPAS_pool_get_array(hexdyn, "solveStress", solveStress)
    call MPAS_pool_get_array(hexdyn, "solveVelocity", solveVelocity)
    call MPAS_pool_get_array(hexdyn, "icePressure", icePressure)
    call MPAS_pool_get_array(hexdyn, "airStressCellU", airStressCellU)
    call MPAS_pool_get_array(hexdyn, "airStressCellV", airStressCellV)
    call MPAS_pool_get_array(hexdyn, "airStressVertexU", airStressVertexU)
    call MPAS_pool_get_array(hexdyn, "airStressVertexV", airStressVertexV)
    call MPAS_pool_get_array(hexdyn, "oceanStressU", oceanStressU)
    call MPAS_pool_get_array(hexdyn, "oceanStressV", oceanStressV)
    call MPAS_pool_get_array(hexdyn, "oceanStressCoeff", oceanStressCoeff)
    call MPAS_pool_get_array(hexdyn, "surfaceTiltForceU", surfaceTiltForceU)
    call MPAS_pool_get_array(hexdyn, "surfaceTiltForceV", surfaceTiltForceV)
    call MPAS_pool_get_array(hexdyn, "totalMassVertexfVertex", totalMassVertexfVertex)
    call MPAS_pool_get_array(hexdyn, "principalStress1", principalStress1)
    call MPAS_pool_get_array(hexdyn, "principalStress2", principalStress2)
    call MPAS_pool_get_array(hexdyn, "replacementPressure", replacementPressure)
    call MPAS_pool_get_array(hexdyn, "uOceanVelocityVertex", uOceanVelocityVertex)
    call MPAS_pool_get_array(hexdyn, "vOceanVelocityVertex", vOceanVelocityVertex)
    call MPAS_pool_get_array(hexdyn, "uVelocity", uVelocity)
    call MPAS_pool_get_array(hexdyn, "vVelocity", vVelocity)

    call MPAS_pool_get_array(ocean_forcing, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(ocean_forcing, "vOceanVelocity", vOceanVelocity)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call MPAS_pool_get_array(atmos_forcing, "uAirVelocity", uAirVelocity)
    call MPAS_pool_get_array(atmos_forcing, "vAirVelocity", vAirVelocity)
    call MPAS_pool_get_array(atmos_forcing, "airDensity", airDensity)

    call MPAS_pool_get_array(mesh, "fVertex", fVertex)

    ! aggregate categories for area and volume
    call aggregate_mass_and_area(mesh, &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         totalMassCell)

    ! interpolate area and mass from cells to vertices
    call interpolate_cell_to_vertex(mesh, &
         iceAreaVertex, &
         iceAreaCell)

    call interpolate_cell_to_vertex(mesh, &
         totalMassVertex, &
         totalMassCell)

    ! calculate computational masks
    call stress_calculation_mask_hex(mesh, &
         solveStress, &
         iceAreaCell, &
         totalMassCell)

    call velocity_calculation_mask_hex(mesh, &
         solveVelocity,   &
         iceAreaVertex,   &
         totalMassVertex, &
         interiorVertex)

    ! calculate the ice strength
    call ice_strength(mesh, &     
         icePressure,   &
         iceAreaCell,   &
         iceVolumeCell, &
         solveStress)

    ! calculate the air stress
    call air_stress(&
         airStressCellU, &
         airStressCellV, &
         uAirVelocity,   & 
         vAirVelocity,   &
         airDensity,     &
         iceAreaCell)

    call interpolate_cell_to_vertex(mesh, &
         airStressVertexU, &
         airStressCellU)

    call interpolate_cell_to_vertex(mesh, &
         airStressVertexV, &
         airStressCellV)

    call interpolate_cell_to_vertex(mesh, &
         uOceanVelocityVertex, &
         uOceanVelocity)

    call interpolate_cell_to_vertex(mesh, &
         vOceanVelocityVertex, &
         vOceanVelocity)

    ! calculate the coriolis force coefficient
    call coriolis_force_coefficient(&
         totalMassVertexfVertex, &
         totalMassVertex, &
         fVertex)

    ! calculate the ocean stress
    call ocean_stress(& 
         oceanStressU,         &
         oceanStressV,         &
         uOceanVelocityVertex, & 
         vOceanVelocityVertex, &
         fVertex,              &
         solveVelocity)

    ! calculate the surface tilt force
    call surface_tilt(&
         surfaceTiltForceU,    &
         surfaceTiltForceV,    &
         uOceanVelocityVertex, & 
         vOceanVelocityVertex, &
         totalMassVertex,      &
         fVertex,              &
         solveVelocity)

    ! subcycle the dynamics
    call subcycle_dynamics_hex(mesh, &
                               icestate,   boundary,  &
                               hexdyn,                &
                               hexweak,    hexvar,    &
                               hexwach,    hexpwl,    &
                               dtDynamics, dtElastic, &
                               configs,    nstep)

    elapsedTime = elapsedTime + dt

    call principal_stresses_driver(&
         hexweak,             &
         hexvar,              &
         replacementPressure, &
         configs,             &
         principalStress1,    &
         principalStress2)

    call ocean_stress_final(&
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

    call plot_vertices_planar(mesh, uVelocity, interiorVertex, "uVelocity.txt")

  end subroutine run_dynamics_hex

  !-------------------------------------------------------------

  subroutine principal_stresses_driver(&
       hexweak,             &
       hexvar,              &
       replacementPressure, &
       configs,             &
       principalStress1,    &
       principalStress2)

    use cice_dynamics_shared, only: &
         principal_stresses         

    type (MPAS_pool_type), pointer, intent(in) :: &
         hexweak, &
         hexvar, &
         configs

    real(kind=RKIND), dimension(:), intent(in) :: &
         replacementPressure
    
    real(kind=RKIND), dimension(:), intent(out) :: &
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

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    call MPAS_pool_get_config(configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)

    ! calculate the principal stresses
    if (trim(config_stress_divergence_scheme) == "weak") then

       call MPAS_pool_get_array(hexweak, "stress11", stress11Weak) 
       call MPAS_pool_get_array(hexweak, "stress22", stress22Weak) 
       call MPAS_pool_get_array(hexweak, "stress12", stress12Weak) 

       call principal_stresses(&
            principalStress1, &
            principalStress2, &
            stress11Weak,     &
            stress22Weak,     &
            stress12Weak,     &
            replacementPressure)

    else if (trim(config_stress_divergence_scheme) == "variational") then

       call MPAS_pool_get_array(hexvar, "stress11", stress11Var) 
       call MPAS_pool_get_array(hexvar, "stress22", stress22Var)
       call MPAS_pool_get_array(hexvar, "stress12", stress12Var)

       call principal_stresses(&
            principalStress1, &
            principalStress2, &
            stress11Var(1,:), &
            stress22Var(1,:), &
            stress12Var(1,:), &
            replacementPressure)

    endif

  end subroutine principal_stresses_driver

  !-------------------------------------------------------------

  subroutine subcycle_dynamics_hex(mesh, &
                                  icestate,   boundary,  &
                                  hexdyn,                &
                                  hexweak,    hexvar,    &
                                  hexwach,    hexpwl,    &
                                  dtDynamics, dtElastic, &
                                  configs,    nstep)

    type(MPAS_pool_type), pointer, intent(in) :: mesh
    type(MPAS_pool_type), pointer, intent(in) :: configs
    type(MPAS_pool_type), pointer :: icestate
    type(MPAS_pool_type), pointer :: boundary
    type(MPAS_pool_type), pointer :: hexdyn
    type(MPAS_pool_type), pointer :: hexweak
    type(MPAS_pool_type), pointer :: hexvar
    type(MPAS_pool_type), pointer :: hexwach
    type(MPAS_pool_type), pointer :: hexpwl

    real(kind=RKIND), intent(in) :: &
         dtDynamics, &
         dtElastic
    
    integer, intent(in) :: &
         nstep

    integer, pointer :: &
         config_elastic_subcycle_number

    integer :: &
         iElasticSubcycle

    call MPAS_pool_get_config(configs, "config_elastic_subcycle_number", config_elastic_subcycle_number)

    do iElasticSubcycle = 1, config_elastic_subcycle_number

       call single_subcycle_dynamics_hex(mesh, &
                                         icestate,   boundary,  &
                                         hexdyn,                &
                                         hexweak,    hexvar,    &
                                         hexwach,    hexpwl,    &
                                         dtDynamics, dtElastic, &
                                         configs,               &
                                         iElasticSubcycle, nstep)

    enddo

  end subroutine subcycle_dynamics_hex

  !-------------------------------------------------------------

  subroutine single_subcycle_dynamics_hex(mesh, &
                                          icestate,   boundary,  &
                                          hexdyn,                &
                                          hexweak,    hexvar,    &
                                          hexwach,    hexpwl,    &
                                          dtDynamics, dtElastic, &
                                          configs,               &
                                          iElasticSubcycle, nstep)

    use cice_mesh, only: interpolate_cell_to_vertex                                    

    use cice_dynamics_shared, only: ocean_stress_coefficient, &
                                    solve_velocity

    use cice_dyn_hex_weak, only: &
         internal_stress_weak

    use cice_dyn_hex_wachspress, only: &
         internal_stress_wachspress

    use cice_dyn_hex_pwl, only: &
         internal_stress_pwl

    use cice_diagnostics, only: &
         plot_cells_spherical, &
         plot_vertices_spherical, &
         plot_pop_grid_vertices, &
         plot_pop_grid_cells

    type(MPAS_pool_type), pointer, intent(in) :: mesh
    type(MPAS_pool_type), pointer, intent(in) :: configs
    type(MPAS_pool_type), pointer :: icestate
    type(MPAS_pool_type), pointer :: boundary
    type(MPAS_pool_type), pointer :: hexdyn
    type(MPAS_pool_type), pointer :: hexweak
    type(MPAS_pool_type), pointer :: hexvar
    type(MPAS_pool_type), pointer :: hexwach
    type(MPAS_pool_type), pointer :: hexpwl

    real(kind=RKIND), intent(in) :: &
         dtDynamics, &
         dtElastic

    integer, intent(in) :: &
         nstep, & !! testing
         iElasticSubcycle !! testing

    integer :: iVertex !! testing

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme, &
         config_variational_basis

    logical, pointer :: &
         config_evp_damping

    integer, pointer :: &
         nVertices

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveStress

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaVertex, &
         totalMassVertex, &
         fVertex, &
         uVelocity, &
         vVelocity, &
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

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(hexdyn, "solveVelocity", solveVelocity)
    call MPAS_pool_get_array(hexdyn, "solveStress", solveStress)
    call MPAS_pool_get_array(hexdyn, "uVelocity", uVelocity)
    call MPAS_pool_get_array(hexdyn, "vVelocity", vVelocity)
    call MPAS_pool_get_array(hexdyn, "icePressure", icePressure)
    call MPAS_pool_get_array(hexdyn, "replacementPressure", replacementPressure)
    call MPAS_pool_get_array(hexdyn, "stressDivergenceU", stressDivergenceU)
    call MPAS_pool_get_array(hexdyn, "stressDivergenceV", stressDivergenceV)
    call MPAS_pool_get_array(hexdyn, "airStressVertexU", airStressVertexU)
    call MPAS_pool_get_array(hexdyn, "airStressVertexV", airStressVertexV)
    call MPAS_pool_get_array(hexdyn, "surfaceTiltForceU", surfaceTiltForceU)
    call MPAS_pool_get_array(hexdyn, "surfaceTiltForceV", surfaceTiltForceV)
    call MPAS_pool_get_array(hexdyn, "totalMassVertexfVertex", totalMassVertexfVertex)
    call MPAS_pool_get_array(hexdyn, "oceanStressU", oceanStressU)
    call MPAS_pool_get_array(hexdyn, "oceanStressV", oceanStressV)
    call MPAS_pool_get_array(hexdyn, "oceanStressCoeff", oceanStressCoeff)
    call MPAS_pool_get_array(hexdyn, "uOceanVelocityVertex", uOceanVelocityVertex)
    call MPAS_pool_get_array(hexdyn, "vOceanVelocityVertex", vOceanVelocityVertex)

    if (trim(config_stress_divergence_scheme) == "weak") then

       call internal_stress_weak(&
            mesh,                hexweak, &
            uVelocity,           vVelocity, &
            icePressure,         replacementPressure, &
            dtElastic,           config_evp_damping, &
            stressDivergenceU,   stressDivergenceV, &
            solveStress,         solveVelocity)

    else if (trim(config_stress_divergence_scheme) == "variational" .and. &
             trim(config_variational_basis) == "wachspress") then

       call internal_stress_wachspress(&
            mesh, &
            hexvar,              hexwach, &
            uVelocity,           vVelocity, &
            icePressure,         replacementPressure, &
            dtElastic,           config_evp_damping, &
            stressDivergenceU,   stressDivergenceV, &
            solveStress,         solveVelocity)

    else if (trim(config_stress_divergence_scheme) == "variational" .and. &
             trim(config_variational_basis) == "pwl") then

       call internal_stress_pwl(&
            mesh, &
            hexvar,              hexpwl, &
            uVelocity,           vVelocity, &
            icePressure,         replacementPressure, &
            dtElastic,           config_evp_damping, &
            stressDivergenceU,   stressDivergenceV, &
            solveStress,         solveVelocity)

    endif

    !stressDivergenceU = 0.0_RKIND
    !stressDivergenceV = 0.0_RKIND

    call ocean_stress_coefficient(& 
         oceanStressCoeff,     &
         uOceanVelocityVertex, & 
         vOceanVelocityVertex, &
         uVelocity,            &   
         vVelocity,            &
         iceAreaVertex,        &
         solveVelocity)

    call solve_velocity(&
         nVertices,              &
         solveVelocity,          &
         uVelocity,              &   
         vVelocity,              &
         totalMassVertex,        & 
         totalMassVertexfVertex, &
         stressDivergenceU,      & 
         stressDivergenceV,      &
         airStressVertexU,       & 
         airStressVertexV,       &
         surfaceTiltForceU,      &
         surfaceTiltForceV,      &
         oceanStressU,           &
         oceanStressV,           &
         oceanStressCoeff,       &
         dtElastic, mesh)

  end subroutine single_subcycle_dynamics_hex

  !-------------------------------------------------------------

  subroutine aggregate_mass_and_area(&
       mesh, &
       iceAreaCategory, &
       iceVolumeCategory, &
       snowVolumeCategory, &
       iceAreaCell, &
       iceVolumeCell, &
       snowVolumeCell, &
       totalMassCell)

    use cice_constants, only: &
         rhoi, &
         rhos

    type (MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory      

    real(kind=RKIND), dimension(:), intent(out) :: &
         totalMassCell, &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell

    integer :: &
         iCell, &
         iCategory

    integer, pointer :: &
         nCells, &
         nCategories

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    do iCell = 1, nCells

       totalMassCell(iCell)  = 0.0_RKIND
       iceAreaCell(iCell)    = 0.0_RKIND
       iceVolumeCell(iCell)  = 0.0_RKIND
       snowVolumeCell(iCell) = 0.0_RKIND

       do iCategory = 1, nCategories

          iceAreaCell(iCell)    = iceAreaCell(iCell) + iceAreaCategory(1,iCategory,iCell)
          iceVolumeCell(iCell)  = iceVolumeCell(iCell) + iceVolumeCategory(1,iCategory,iCell)
          snowVolumeCell(iCell) = snowVolumeCell(iCell) + snowVolumeCategory(1,iCategory,iCell)

          totalMassCell(iCell)  = totalMassCell(iCell) + &
               iceVolumeCell(iCell) * rhoi + &
               snowVolumeCell(iCell) * rhos

       enddo ! iCategory

    enddo ! iCell

  end subroutine aggregate_mass_and_area

  !-------------------------------------------------------------

  subroutine stress_calculation_mask_hex(mesh,     solveStress, &
                                         areaCell, massCell)

    use cice_dynamics_shared, only: areaMinimum, massMinimum

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    integer, dimension(:), intent(out) :: &
         solveStress

    real(kind=RKIND), dimension(:), intent(in) :: &
         areaCell, &
         massCell

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

    do iCell = 1, nCells

       solveStress(iCell) = 0

       if (areaCell(iCell) > areaMinimum .and. &
           massCell(iCell) > massMinimum) then

          ! this cell has sufficient ice
          solveStress(iCell) = 1

       else

          ! test neighbouring cells to see if have sufficient ice
          do iCellOnCell = 1, nEdgesOnCell(iCell)

             iCellNeighbour = cellsOnCell(iCellOnCell,iCell)

             if (areaCell(iCellNeighbour) > areaMinimum .and. &
                 massCell(iCellNeighbour) > massMinimum) then

                solveStress(iCell) = 1
                exit

             endif

          enddo ! iCellOnCell

       endif

    enddo ! iCell

  end subroutine stress_calculation_mask_hex

  !-------------------------------------------------------------

  subroutine velocity_calculation_mask_hex(mesh,       solveVelocity, &
                                           areaVertex, massVertex, &
                                           interiorVertex)

    use cice_dynamics_shared, only: areaMinimum, massMinimum

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    integer, dimension(:), intent(out) :: &
         solveVelocity

    real(kind=RKIND), dimension(:), intent(in) :: &
         areaVertex, &
         massVertex

    integer, dimension(:), intent(in) :: &
         interiorVertex    

    integer :: &
         iVertex

    integer, pointer :: &
         nVertices

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    do iVertex = 1, nVertices

       solveVelocity(iVertex) = 0

       if (interiorVertex(iVertex) == 1 .and. &
           areaVertex(iVertex) > areaMinimum .and. &
           massVertex(iVertex) > massMinimum) then

          ! this vertex has sufficient ice
          solveVelocity(iVertex) = 1

       endif

    enddo ! iVertex

  end subroutine velocity_calculation_mask_hex

  !-------------------------------------------------------------

end module cice_dynamics_hex

