module cice_dynamics_hex

  use mpas_grid_types

  implicit none

  private
  public :: init_dynamics_hex, &
            run_dynamics_hex
  
contains

  !-------------------------------------------------------------
  ! Initialization
  !-------------------------------------------------------------
  
  subroutine init_dynamics_hex(mesh,     dt,          &
                               boundary, normal,      &
                               hexvar,   hexpwl,      &
                               divergenceFormulation, &
                               gridCurvature)

    use cice_dynamics_shared, only: init_evp
    use cice_dyn_hex_weak, only: init_dynamics_hex_weak
    use cice_dyn_hex_variational, only: init_dynamics_hex_variational
    use cice_dyn_hex_pwl, only: init_dynamics_hex_pwl

    use cice_testing, only: plot_pop_grid_vertices

    type(mesh_type), intent(inout) :: mesh

    type(normal_type),   pointer :: normal
    type(hexvar_type),   pointer :: hexvar
    type(hexpwl_type),   pointer :: hexpwl
    type(boundary_type), pointer :: boundary

    real(kind=RKIND), intent(in) :: dt

    character(len=*), intent(in) :: &
         divergenceFormulation, &
         gridCurvature

    call init_evp(dt)

    if (trim(divergenceFormulation) == "weak") then

       call init_dynamics_hex_weak(mesh, boundary, normal, gridCurvature)

    else if (trim(divergenceFormulation) == "variational") then

       call init_dynamics_hex_variational(mesh, hexvar, boundary)

    else if (trim(divergenceFormulation) == "pwl") then

       call init_dynamics_hex_pwl(mesh, hexvar, hexpwl, boundary)

    endif

    call plot_pop_grid_vertices(mesh, mesh % fVertex % array, boundary % interiorVertex % array, "coriolis.txt", 0)
    call plot_pop_grid_vertices(mesh, mesh % latVertex % array, boundary % interiorVertex % array, "latVertex.txt", 0)

  end subroutine init_dynamics_hex

  !-------------------------------------------------------------
  ! Time stepping
  !-------------------------------------------------------------

  subroutine run_dynamics_hex(mesh, &
                              icestate, boundary,        &
                              normal,   atmos_forcing,   &
                              hexdyn,   hexweak,         &
                              hexvar,   hexpwl,          &
                              hexfor,                    &
                              dt, divergenceFormulation, &
                              evpDamping, &
                              nstep)

    use cice_dynamics_shared, only: ice_strength, &
                                    interpolate_cell_to_vertex, &
                                    air_stress, &
                                    ocean_stress, &
                                    surface_tilt, &
                                    principal_stresses, &
                                    dtDynamics, &
                                    dtElastic, &
                                    elapsedTime

    use cice_testing, only: gnuplot_triangle, gnuplot_cell, writeout_minmax, plot_pop_grid_cells, &
         plot_cells_spherical, &
         plot_vertices_spherical, &
         plot_pop_grid_vertices

    type(mesh_type), intent(in) :: mesh

    type(icestate_type),      pointer :: icestate
    type(boundary_type),      pointer :: boundary
    type(normal_type),        pointer :: normal
    type(atmos_forcing_type), pointer :: atmos_forcing    
    type(hexdyn_type),        pointer :: hexdyn
    type(hexweak_type),       pointer :: hexweak
    type(hexvar_type),        pointer :: hexvar
    type(hexpwl_type),        pointer :: hexpwl
    type(hexfor_type),        pointer :: hexfor

    real(kind=RKIND), intent(in) :: &
         dt

    character(len=*), intent(in) :: &
         divergenceFormulation

    logical, intent(in) :: &
         evpDamping

    integer, intent(in) :: &
         nstep

    ! interpolate area and mass from cells to vertices
    call interpolate_cell_to_vertex(mesh, &
         icestate % iceAreaVertex % array, &
         icestate % iceAreaCell % array)

    call interpolate_cell_to_vertex(mesh, &
         icestate % totalMassVertex % array, &
         icestate % totalMassCell % array)

    call plot_pop_grid_vertices(mesh, icestate % totalMassVertex % array, boundary % interiorVertex % array, "totalMassVertex.txt", nstep)


    ! calculate computational masks
    call stress_calculation_mask_hex(mesh, &
         hexdyn % solveStress % array, &
         icestate % iceAreaCell % array, &
         icestate % totalMassCell % array)

    call velocity_calculation_mask_hex(mesh, &
         hexdyn % solveVelocity % array, &
         icestate % iceAreaVertex % array, &
         icestate % totalMassVertex % array, &
         boundary % interiorVertex % array)



    ! calculate the ice strength
    call ice_strength(mesh, &     
         hexdyn % icePressure % array,     &
         icestate % iceAreaCell % array,   &
         icestate % iceVolumeCell % array, &
         hexdyn % solveStress % array)

    call plot_pop_grid_cells(mesh, hexdyn % icePressure % array, "strength.txt", nstep)

    ! calculate the air stress
    call air_stress(&
         hexdyn % airStressCellU % array,      &
         hexdyn % airStressCellV % array,      &
         atmos_forcing % uAirVelocity % array, & 
         atmos_forcing % vAirVelocity % array, &
         atmos_forcing % airDensity % array,   &
         icestate % iceAreaCell % array)

    call plot_pop_grid_cells(mesh, atmos_forcing % uAirVelocity % array, "uAirVelocity.txt", nstep)
    call plot_pop_grid_cells(mesh, atmos_forcing % vAirVelocity % array, "vAirVelocity.txt", nstep)
    call plot_pop_grid_cells(mesh, atmos_forcing % uAirVelocityPOP % array, "uAirVelocityPOP.txt", nstep)
    call plot_pop_grid_cells(mesh, atmos_forcing % vAirVelocityPOP % array, "vAirVelocityPOP.txt", nstep)
    call plot_pop_grid_cells(mesh, atmos_forcing % airTemperature % array, "airTemperature.txt", nstep)
    call plot_pop_grid_cells(mesh, atmos_forcing % airDensity % array, "airDensity.txt", nstep)
    call plot_pop_grid_cells(mesh, hexdyn % airStressCellU % array, "airStressCellU.txt", nstep)
    call plot_pop_grid_cells(mesh, hexdyn % airStressCellV % array, "airStressCellV.txt", nstep)


    call interpolate_cell_to_vertex(mesh, &
         hexdyn % airStressVertexU % array, &
         hexdyn % airStressCellU % array)

    call interpolate_cell_to_vertex(mesh, &
         hexdyn % airStressVertexV % array, &
         hexdyn % airStressCellV % array)

    call plot_pop_grid_vertices(mesh, hexdyn % airStressVertexU % array, boundary % interiorVertex % array, "airStressVertexU.txt", nstep)
    call plot_pop_grid_vertices(mesh, hexdyn % airStressVertexV % array, boundary % interiorVertex % array, "airStressVertexV.txt", nstep)

    ! calculate the ocean stress
    call ocean_stress(& 
         hexdyn % oceanStressU % array,   &
         hexdyn % oceanStressV % array,   &
         hexfor % uOceanVelocity % array, & 
         hexfor % vOceanVelocity % array, &
         mesh % fVertex % array,          &
         hexdyn % solveVelocity % array)

    !call plot_pop_grid_vertices(mesh, hexfor % uOceanVelocity % array, boundary % interiorVertex % array, "uOceanVelocity.txt", nstep)
    !call plot_pop_grid_vertices(mesh, hexfor % uOceanVelocity % array, boundary % interiorVertex % array, "uOceanVelocity.txt", nstep)
    !call plot_pop_grid_vertices(mesh, hexdyn % oceanStressU % array, boundary % interiorVertex % array, "oceanStressU.txt", nstep)
    !call plot_pop_grid_vertices(mesh, hexdyn % oceanStressV % array, boundary % interiorVertex % array, "oceanStressV.txt", nstep)

    ! calculate the surface tilt force
    call surface_tilt(&
         hexdyn % surfaceTiltForceU % array, &
         hexdyn % surfaceTiltForceV % array, &
         hexfor % uOceanVelocity % array,    & 
         hexfor % vOceanVelocity % array,    &
         icestate % totalMassVertex % array, &
         mesh % fVertex % array,             &
         hexdyn % solveVelocity % array)

    !call plot_pop_grid_vertices(mesh, hexdyn % surfaceTiltForceU % array, boundary % interiorVertex % array, "surfaceTiltForceU.txt", nstep)
    !call plot_pop_grid_vertices(mesh, hexdyn % surfaceTiltForceV % array, boundary % interiorVertex % array, "surfaceTiltForceV.txt", nstep)

    ! subcycle the dynamics
    call subcycle_dynamics_hex(mesh, &
                               icestate,   boundary,  &
                               normal,     hexdyn,    &
                               hexweak,    hexvar,    &
                               hexpwl,     hexfor,    &
                               dtDynamics, dtElastic, &
                               divergenceFormulation, &
                               evpDamping, nstep)

    elapsedTime = elapsedTime + dt

    ! calculate the principal stresses
    if (trim(divergenceFormulation) == "weak") then

       call principal_stresses(&
            hexdyn % principalStress1 % array, &
            hexdyn % principalStress2 % array, &
            hexweak % stress11 % array, &
            hexweak % stress22 % array, &
            hexweak % stress12 % array, &
            hexdyn % replacementPressure % array)

    else if (trim(divergenceFormulation) == "variational") then

       call principal_stresses(&
            hexdyn % principalStress1 % array, &
            hexdyn % principalStress2 % array, &
            hexvar % stress11 % array(1,:), &
            hexvar % stress22 % array(1,:), &
            hexvar % stress12 % array(1,:), &
            hexdyn % replacementPressure % array)

    endif

  end subroutine run_dynamics_hex

  !-------------------------------------------------------------

  subroutine subcycle_dynamics_hex(mesh, &
                                  icestate,   boundary,  &
                                  normal,     hexdyn,    &
                                  hexweak,    hexvar,    &
                                  hexpwl,     hexfor,    &
                                  dtDynamics, dtElastic, &
                                  divergenceFormulation, &
                                  evpDamping, nstep)

    use cice_dynamics_shared, only: nElasticSubcycle

    type(mesh_type), intent(in) :: mesh

    type(icestate_type), pointer :: icestate
    type(boundary_type), pointer :: boundary
    type(normal_type),   pointer :: normal
    type(hexdyn_type),   pointer :: hexdyn
    type(hexweak_type),  pointer :: hexweak
    type(hexvar_type),   pointer :: hexvar
    type(hexpwl_type),   pointer :: hexpwl
    type(hexfor_type),   pointer :: hexfor

    real(kind=RKIND), intent(in) :: &
         dtDynamics, &
         dtElastic

    character(len=*), intent(in) :: &
         divergenceFormulation

    logical, intent(in) :: &
         evpDamping
    
    integer, intent(in) :: &
         nstep

    integer :: &
         iElasticSubcycle

    do iElasticSubcycle = 1, 10!nElasticSubcycle

       write(*,*) "subcycle: ", iElasticSubcycle
       call single_subcycle_dynamics_hex(mesh, &
                                         icestate,   boundary,  &
                                         normal,     hexdyn,    &
                                         hexweak,    hexvar,    &
                                         hexpwl,     hexfor,    &
                                         dtDynamics, dtElastic, &
                                         divergenceFormulation, &
                                         evpDamping,            &
                                         iElasticSubcycle, nstep)

    enddo

    stop

  end subroutine subcycle_dynamics_hex

  !-------------------------------------------------------------

  subroutine single_subcycle_dynamics_hex(mesh, &
                                          icestate,   boundary,  &
                                          normal,     hexdyn,    &
                                          hexweak,    hexvar,    &
                                          hexpwl,     hexfor,    &
                                          dtDynamics, dtElastic, &
                                          divergenceFormulation, &
                                          evpDamping,            &
                                          iElasticSubcycle, nstep)

    use cice_dynamics_shared, only: interpolate_cell_to_vertex, &
                                    ocean_stress_coefficient, &
                                    solve_velocity

    use cice_dyn_hex_weak, only: strain_tensor_hex_weak, &
                                 stress_tensor_hex_weak, &
                                 stress_divergence_hex_weak, &
                                 stress_divergence_hex_weak_corners

    use cice_dyn_hex_variational, only: strain_tensor_hex_variational, &
                                        stress_tensor_hex_variational, &
                                        stress_divergence_hex_variational

    use cice_dyn_hex_pwl, only: strain_tensor_hex_pwl, &
                                stress_divergence_hex_pwl

    use cice_testing, only: &
         plot_cells_spherical, &
         plot_vertices_spherical, &
         plot_pop_grid_vertices, &
         plot_pop_grid_cells

    type(mesh_type), intent(in) :: mesh

    type(icestate_type), pointer :: icestate
    type(boundary_type), pointer :: boundary
    type(normal_type),   pointer :: normal
    type(hexdyn_type),   pointer :: hexdyn
    type(hexweak_type),  pointer :: hexweak
    type(hexvar_type),   pointer :: hexvar
    type(hexpwl_type),   pointer :: hexpwl
    type(hexfor_type),   pointer :: hexfor

    real(kind=RKIND), intent(in) :: &
         dtDynamics, &
         dtElastic

    logical, intent(in) :: &
         evpDamping

    integer, intent(in) :: &
         nstep, & !! testing
         iElasticSubcycle !! testing

    integer :: iVertex !! testing

    character(len=*), intent(in) :: divergenceFormulation

    do iVertex = 1, mesh % nVertices
    
       if (boundary % interiorVertex % array(iVertex) == 1) then
       
          hexdyn % uVelocity % array(iVertex) = sin(4.0_RKIND * mesh % latVertex % array(iVertex) + 2.0_RKIND * 3.142_RKIND * sin(2.0_RKIND * mesh % lonVertex % array(iVertex)))
          hexdyn % vVelocity % array(iVertex) = 0.0_RKIND
          
       else

          hexdyn % uVelocity % array(iVertex) = 0.0_RKIND
          hexdyn % vVelocity % array(iVertex) = 0.0_RKIND

       endif

    enddo ! iVertex

    if (trim(divergenceFormulation) == "weak") then

       call strain_tensor_hex_weak(mesh, &
            hexweak % strain11 % array,           &
            hexweak % strain22 % array,           &
            hexweak % strain12 % array,           &
            hexdyn % uVelocity % array,           &
            hexdyn % vVelocity % array,           &
            normal % normalVectorPolygon % array, &
            normal % latCellRotated % array,      &
            hexdyn % solveStress % array)

       call plot_pop_grid_cells(mesh, hexweak % strain11 % array, "strain11.txt", iElasticSubcycle)
       call plot_pop_grid_cells(mesh, hexweak % strain22 % array, "strain22.txt", iElasticSubcycle)
       call plot_pop_grid_cells(mesh, hexweak % strain12 % array, "strain12.txt", iElasticSubcycle)

       call stress_tensor_hex_weak(mesh, &
            hexweak % stress11 % array,           &
            hexweak % stress22 % array,           &
            hexweak % stress12 % array,           &
            hexweak % strain11 % array,           &
            hexweak % strain22 % array,           &
            hexweak % strain12 % array,           &
            hexdyn % icePressure % array,         &
            hexdyn % replacementPressure % array, &
            hexdyn % solveStress % array,          &       
            dtElastic, evpDamping)

       call stress_divergence_hex_weak(mesh, &
            hexdyn % stressDivergenceU % array,    &
            hexdyn % stressDivergenceV % array,    &
            hexweak % stress11 % array,            &
            hexweak % stress22 % array,            &
            hexweak % stress12 % array,            &
            normal % normalVectorTriangle % array, &
            normal % latVertexRotated % array,     &
            hexdyn % solveVelocity % array)

       call plot_pop_grid_vertices(mesh, hexdyn % stressDivergenceU % array, boundary % interiorVertex % array, "divstressu.txt", iElasticSubcycle)
       call plot_pop_grid_vertices(mesh, hexdyn % stressDivergenceV % array, boundary % interiorVertex % array, "divstressv.txt", iElasticSubcycle)

    else if (trim(divergenceFormulation) == "variational") then

       call strain_tensor_hex_variational(mesh, &
            hexvar % strain11 % array,       &
            hexvar % strain22 % array,       &
            hexvar % strain12 % array,       &
            hexdyn % uVelocity % array,      &
            hexdyn % vVelocity % array,      &
            hexvar % basisGradientU % array, &
            hexvar % basisGradientV % array, &
            hexdyn % solveStress % array)

       call stress_tensor_hex_variational(mesh, &
            hexvar % stress11 % array,            &
            hexvar % stress22 % array,            &
            hexvar % stress12 % array,            &
            hexvar % strain11 % array,            &
            hexvar % strain22 % array,            &
            hexvar % strain12 % array,            &
            hexdyn % icePressure % array,         &
            hexdyn % replacementPressure % array, &
            hexdyn % solveStress % array,          &
            dtElastic, evpDamping)

       call stress_divergence_hex_variational(mesh, &
            hexdyn % stressDivergenceU % array,    &
            hexdyn % stressDivergenceV % array,    &
            hexvar % stress11 % array,             &
            hexvar % stress22 % array,             &
            hexvar % stress12 % array,             &
            hexvar % basisIntegralsU % array,      &
            hexvar % basisIntegralsV % array,      &
            hexvar % cellVerticesAtVertex % array, &
            hexdyn % solveVelocity % array)

    else if (trim(divergenceFormulation) == "pwl") then

       call strain_tensor_hex_pwl(mesh, &
            hexvar % strain11 % array,       &
            hexvar % strain22 % array,       &
            hexvar % strain12 % array,       &
            hexdyn % uVelocity % array,      &
            hexdyn % vVelocity % array,      &
            hexpwl % basisGradientU % array, &
            hexpwl % basisGradientV % array, &
            hexdyn % solveStress % array)

       call stress_tensor_hex_variational(mesh, &
            hexvar % stress11 % array,            &
            hexvar % stress22 % array,            &
            hexvar % stress12 % array,            &
            hexvar % strain11 % array,            &
            hexvar % strain22 % array,            &
            hexvar % strain12 % array,            &
            hexdyn % icePressure % array,         &
            hexdyn % replacementPressure % array, &
            hexdyn % solveStress % array,          &
            dtElastic, evpDamping)

       call stress_divergence_hex_pwl(mesh, &
            hexdyn % stressDivergenceU % array,    &
            hexdyn % stressDivergenceV % array,    &
            hexvar % stress11 % array,             &
            hexvar % stress22 % array,             &
            hexvar % stress12 % array,             &
            hexpwl % basisGradientU % array,       &
            hexpwl % basisGradientV % array,       &
            hexpwl % basisSubArea % array,         &
            hexvar % cellVerticesAtVertex % array, &
            hexdyn % solveVelocity % array)

    endif

    call ocean_stress_coefficient(& 
         hexdyn % oceanStressCoeff % array, &
         hexfor % uOceanVelocity % array,   & 
         hexfor % vOceanVelocity % array,   &
         hexdyn % uVelocity % array,        &   
         hexdyn % vVelocity % array,        &
         icestate % iceAreaVertex % array,  &
         hexdyn % solveVelocity % array)

    !call solve_velocity(mesh % nVertices, &
    !     hexdyn % solveVelocity % array,      &
    !     hexdyn % uVelocity % array,         &   
    !     hexdyn % vVelocity % array,         &
    !     icestate % totalMassVertex % array, & 
    !     mesh % fVertex % array,             &
    !     hexdyn % stressDivergenceU % array, & 
    !     hexdyn % stressDivergenceV % array, &
    !     hexdyn % airStressVertexU % array,  & 
    !     hexdyn % airStressVertexV % array,  &
    !     hexdyn % surfaceTiltForceU % array, &
    !     hexdyn % surfaceTiltForceV % array, &
    !     hexdyn % oceanStressU % array,      &
    !     hexdyn % oceanStressV % array,      &
    !     hexdyn % oceanStressCoeff % array,  &
    !     dtElastic,mesh)

    call plot_pop_grid_vertices(mesh, hexdyn % uVelocity % array, boundary % interiorVertex % array, "uvel.txt", iElasticSubcycle)
    call plot_pop_grid_vertices(mesh, hexdyn % vVelocity % array, boundary % interiorVertex % array, "vvel.txt", iElasticSubcycle)

  end subroutine single_subcycle_dynamics_hex

  !-------------------------------------------------------------

  subroutine stress_calculation_mask_hex(mesh,     solveStress, &
                                         areaCell, massCell)

    use cice_dynamics_shared, only: areaMinimum, massMinimum

    type(mesh_type), intent(in) :: mesh

    integer, dimension(:), intent(out) :: &
         solveStress

    real(kind=RKIND), dimension(:), intent(in) :: &
         areaCell, &
         massCell

    integer :: &
         iCell, &
         iCellOnCell, &
         iCellNeighbour

    do iCell = 1, mesh % nCells

       solveStress(iCell) = 0

       if (areaCell(iCell) > areaMinimum .and. &
           massCell(iCell) > massMinimum) then

          ! this cell has sufficient ice
          solveStress(iCell) = 1

       else

          ! test neighbouring cells to see if have sufficient ice
          do iCellOnCell = 1, mesh % nEdgesOnCell % array(iCell)

             iCellNeighbour = mesh % cellsOnCell % array(iCellOnCell,iCell)

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

    type(mesh_type), intent(in) :: mesh

    integer, dimension(:), intent(out) :: &
         solveVelocity

    real(kind=RKIND), dimension(:), intent(in) :: &
         areaVertex, &
         massVertex

    integer, dimension(:), intent(in) :: &
         interiorVertex    

    integer :: &
         iVertex

    do iVertex = 1, mesh % nVertices

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

