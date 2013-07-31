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
                              evpDamping)

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

    ! interpolate area and mass from cells to vertices
    call interpolate_cell_to_vertex(mesh, &
         icestate % iceAreaVertex % array, &
         icestate % iceAreaCell % array)

    call interpolate_cell_to_vertex(mesh, &
         icestate % totalMassVertex % array, &
         icestate % totalMassCell % array)
 
    call plot_cells_spherical(mesh, &
         icestate % totalMassCell % array, &
         "totalMassCell.txt")

    !call plot_vertices_spherical(mesh, &
    !     icestate % totalMassVertex % array, &
    !     boundary % interiorVertex % array, &
    !     "totalMassVertex.txt")


    call plot_pop_grid_vertices(mesh, &
         icestate % totalMassVertex % array, &
         boundary % interiorVertex % array, &
         "totalMassVertex.txt")


    stop

    ! calculate computational masks
    call stress_calculation_mask_hex(mesh, &
         hexdyn % solveStress % array, &
         icestate % iceAreaCell % array, &
         icestate % totalMassCell % array)

    call plot_cells_spherical(mesh, &
         real(hexdyn % solveStress % array,RKIND), &
         "solveStress.txt")

    call velocity_calculation_mask_hex(mesh, &
         hexdyn % solveVelocity % array, &
         icestate % iceAreaVertex % array, &
         icestate % totalMassVertex % array, &
         boundary % interiorVertex % array)

    call plot_vertices_spherical(mesh, &
         real(hexdyn % solveVelocity % array,RKIND), &
         boundary % interiorVertex % array, &
         "solveVelocity.txt")

    ! calculate the ice strength
    call ice_strength(mesh, &     
         hexdyn % icePressure % array,     &
         icestate % iceAreaCell % array,   &
         icestate % iceVolumeCell % array, &
         hexdyn % solveStress % array)

    call plot_cells_spherical(mesh, &
         real(hexdyn % icePressure % array,RKIND), &
         "icePressure.txt")

    ! calculate the air stress
    call air_stress(&
         hexdyn % airStressCellU % array,      &
         hexdyn % airStressCellV % array,      &
         atmos_forcing % uAirVelocity % array, & 
         atmos_forcing % vAirVelocity % array, &
         icestate % iceAreaCell % array)

    call plot_cells_spherical(mesh, &
         real(hexdyn % airStressCellU % array,RKIND), &
         "airStressCellU.txt")

    call plot_cells_spherical(mesh, &
         real(hexdyn % airStressCellV % array,RKIND), &
         "airStressCellV.txt")
  
    call interpolate_cell_to_vertex(mesh, &
         hexdyn % airStressVertexU % array, &
         hexdyn % airStressCellU % array)

    call plot_vertices_spherical(mesh, &
         hexdyn % airStressVertexU % array, &
         boundary % interiorVertex % array, &
         "airStressVertexU.txt")

    call interpolate_cell_to_vertex(mesh, &
         hexdyn % airStressVertexV % array, &
         hexdyn % airStressCellV % array)

    call plot_vertices_spherical(mesh, &
         hexdyn % airStressVertexV % array, &
         boundary % interiorVertex % array, &
         "airStressVertexV.txt")

    ! calculate the ocean stress
    call ocean_stress(& 
         hexdyn % oceanStressU % array,   &
         hexdyn % oceanStressV % array,   &
         hexfor % uOceanVelocity % array, & 
         hexfor % vOceanVelocity % array, &
         mesh % fVertex % array,          &
         hexdyn % solveVelocity % array)

    call plot_vertices_spherical(mesh, &
         hexdyn % oceanStressU % array, &
         boundary % interiorVertex % array, &
         "oceanStressU.txt")

    call plot_vertices_spherical(mesh, &
         hexdyn % oceanStressV % array, &
         boundary % interiorVertex % array, &
         "oceanStressV.txt")

    ! calculate the surface tilt force
    call surface_tilt(&
         hexdyn % surfaceTiltForceU % array, &
         hexdyn % surfaceTiltForceV % array, &
         hexfor % uOceanVelocity % array,    & 
         hexfor % vOceanVelocity % array,    &
         icestate % totalMassVertex % array, &
         mesh % fVertex % array,             &
         hexdyn % solveVelocity % array)

    call plot_vertices_spherical(mesh, &
         hexdyn % surfaceTiltForceU % array, &
         boundary % interiorVertex % array, &
         "surfaceTiltForceU.txt")

    call plot_vertices_spherical(mesh, &
         hexdyn % surfaceTiltForceV % array, &
         boundary % interiorVertex % array, &
         "surfaceTiltForceV.txt")

    ! subcycle the dynamics
    call subcycle_dynamics_hex(mesh, &
                               icestate,   boundary,  &
                               normal,     hexdyn,    &
                               hexweak,    hexvar,    &
                               hexpwl,     hexfor,    &
                               dtDynamics, dtElastic, &
                               divergenceFormulation, &
                               evpDamping)

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


    call plot_vertices_spherical(mesh, &
         hexdyn % uVelocity % array, &
         boundary % interiorVertex % array, &
         "uVelocity.txt")
    
    call plot_vertices_spherical(mesh, &
         hexdyn % vVelocity % array, &
         boundary % interiorVertex % array, &
         "vVelocity.txt")


    stop

  end subroutine run_dynamics_hex

  !-------------------------------------------------------------

  subroutine subcycle_dynamics_hex(mesh, &
                                  icestate,   boundary,  &
                                  normal,     hexdyn,    &
                                  hexweak,    hexvar,    &
                                  hexpwl,     hexfor,    &
                                  dtDynamics, dtElastic, &
                                  divergenceFormulation, &
                                  evpDamping)

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

    integer :: &
         iElasticSubcycle

    do iElasticSubcycle = 1, nElasticSubcycle

       write(*,*) "subcycle: ", iElasticSubcycle
       call single_subcycle_dynamics_hex(mesh, &
                                         icestate,   boundary,  &
                                         normal,     hexdyn,    &
                                         hexweak,    hexvar,    &
                                         hexpwl,     hexfor,    &
                                         dtDynamics, dtElastic, &
                                         divergenceFormulation, &
                                         evpDamping)

    enddo

  end subroutine subcycle_dynamics_hex

  !-------------------------------------------------------------

  subroutine single_subcycle_dynamics_hex(mesh, &
                                          icestate,   boundary,  &
                                          normal,     hexdyn,    &
                                          hexweak,    hexvar,    &
                                          hexpwl,     hexfor,    &
                                          dtDynamics, dtElastic, &
                                          divergenceFormulation, &
                                          evpDamping)

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
         plot_vertices_spherical

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

    character(len=*), intent(in) :: divergenceFormulation

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

       !call plot_cells_spherical(mesh, &
       !     hexweak % strain11 % array, &
       !     "strain11.txt")

       !call plot_cells_spherical(mesh, &
       !     hexweak % strain22 % array, &
       !     "strain22.txt")

       !call plot_cells_spherical(mesh, &
       !     hexweak % strain12 % array, &
       !     "strain12.txt")



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

       !call plot_cells_spherical(mesh, &
       !     hexweak % stress11 % array, &
       !     "stress11.txt")

       !call plot_cells_spherical(mesh, &
       !     hexweak % stress22 % array, &
       !     "stress22.txt")

       !call plot_cells_spherical(mesh, &
       !     hexweak % stress12 % array, &
       !     "stress12.txt")




       call stress_divergence_hex_weak(mesh, &
            hexdyn % stressDivergenceU % array,    &
            hexdyn % stressDivergenceV % array,    &
            hexweak % stress11 % array,            &
            hexweak % stress22 % array,            &
            hexweak % stress12 % array,            &
            normal % normalVectorTriangle % array, &
            normal % latVertexRotated % array,     &
            hexdyn % solveVelocity % array)


       !call plot_vertices_spherical(mesh, &
       !     hexdyn % stressDivergenceU % array, &
       !     boundary % interiorVertex % array, &
       !     "stressDivergenceU.txt")
       
       !call plot_vertices_spherical(mesh, &
       !     hexdyn % stressDivergenceV % array, &
       !     boundary % interiorVertex % array, &
       !     "stressDivergenceV.txt")





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

    !call plot_vertices_spherical(mesh, &
    !     hexdyn % oceanStressCoeff % array, &
    !     boundary % interiorVertex % array, &
    !     "oceanStressCoeff.txt")
    



    call solve_velocity(mesh % nVertices, &
         hexdyn % solveVelocity % array,      &
         hexdyn % uVelocity % array,         &   
         hexdyn % vVelocity % array,         &
         icestate % totalMassVertex % array, & 
         mesh % fVertex % array,             &
         hexdyn % stressDivergenceU % array, & 
         hexdyn % stressDivergenceV % array, &
         hexdyn % airStressVertexU % array,  & 
         hexdyn % airStressVertexV % array,  &
         hexdyn % surfaceTiltForceU % array, &
         hexdyn % surfaceTiltForceV % array, &
         hexdyn % oceanStressU % array,      &
         hexdyn % oceanStressV % array,      &
         hexdyn % oceanStressCoeff % array,  &
         dtElastic)

    !call plot_vertices_spherical(mesh, &
    !     hexdyn % uVelocity % array, &
    !     boundary % interiorVertex % array, &
    !     "uVelocity.txt")
    
    !call plot_vertices_spherical(mesh, &
    !     hexdyn % vVelocity % array, &
    !     boundary % interiorVertex % array, &
    !     "vVelocity.txt")


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

