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
  
  subroutine init_dynamics_hex(mesh,     dt,          &
                               boundary, normal,      &
                               hexvar,                &
                               hexwach,   hexpwl,     &
                               divergenceFormulation, &
                               gridCurvature)

    use cice_dynamics_shared, only: init_evp
    use cice_dyn_hex_weak, only: init_dynamics_hex_weak
    use cice_dyn_hex_wachspress, only: init_dynamics_hex_wachspress
    use cice_dyn_hex_pwl, only: init_dynamics_hex_pwl

    use cice_diagnostics, only: plot_pop_grid_vertices

    type(MPAS_pool_type), pointer, intent(inout) :: mesh

    type(MPAS_pool_type), pointer :: normal
    type(MPAS_pool_type), pointer :: hexvar
    type(MPAS_pool_type), pointer :: hexwach
    type(MPAS_pool_type), pointer :: hexpwl
    type(MPAS_pool_type), pointer :: boundary

    real(kind=RKIND), intent(in) :: dt

    character(len=*), intent(in) :: &
         divergenceFormulation, &
         gridCurvature
    
    call init_evp(dt)

    if (trim(divergenceFormulation) == "weak") then

       call init_dynamics_hex_weak(mesh, boundary, normal, gridCurvature)

    else if (trim(divergenceFormulation) == "variational") then

       call init_dynamics_hex_wachspress(mesh, hexvar, hexwach, boundary)

    else if (trim(divergenceFormulation) == "pwl") then

       call init_dynamics_hex_pwl(mesh, hexvar, hexpwl, boundary)

    endif

  end subroutine init_dynamics_hex

  !-------------------------------------------------------------
  ! Time stepping
  !-------------------------------------------------------------

  subroutine run_dynamics_hex(mesh, &
                              tracers,  icestate,        &
                              boundary,                  &
                              normal,   atmos_forcing,   &
                              hexdyn,   hexweak,         &
                              hexvar,                    &
                              hexwach,  hexpwl,          &
                              hexfor,                    &
                              dt, divergenceFormulation, &
                              evpDamping, &
                              nstep)

    use cice_mesh, only: interpolate_cell_to_vertex

    use cice_dynamics_shared, only: ice_strength, &           
                                    air_stress, &
                                    ocean_stress, &
                                    surface_tilt, &
                                    dtDynamics, &
                                    dtElastic, &
                                    elapsedTime

    use cice_diagnostics, only: gnuplot_triangle, gnuplot_cell, writeout_minmax, plot_pop_grid_cells, &
         plot_cells_spherical, &
         plot_vertices_spherical, &
         plot_pop_grid_vertices

    type(MPAS_pool_type), pointer, intent(in) :: mesh
    type(MPAS_pool_type), pointer :: tracers
    type(MPAS_pool_type), pointer :: icestate
    type(MPAS_pool_type), pointer :: boundary
    type(MPAS_pool_type), pointer :: normal
    type(MPAS_pool_type), pointer :: atmos_forcing    
    type(MPAS_pool_type), pointer :: hexdyn
    type(MPAS_pool_type), pointer :: hexweak
    type(MPAS_pool_type), pointer :: hexvar
    type(MPAS_pool_type), pointer :: hexwach
    type(MPAS_pool_type), pointer :: hexpwl
    type(MPAS_pool_type), pointer :: hexfor

    real(kind=RKIND), intent(in) :: &
         dt

    character(len=*), intent(in) :: &
         divergenceFormulation

    logical, intent(in) :: &
         evpDamping

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
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         principalStress1, &
         principalStress2, &
         replacementPressure, &
         uOceanVelocity, &
         vOceanVelocity, &
         airDensity, &
         uAirVelocity, &
         vAirVelocity, &
         fVertex

    real(kind=RKIND), dimension(:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    ! init variables
    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory)

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
    call MPAS_pool_get_array(hexdyn, "surfaceTiltForceU", surfaceTiltForceU)
    call MPAS_pool_get_array(hexdyn, "surfaceTiltForceV", surfaceTiltForceV)
    call MPAS_pool_get_array(hexdyn, "principalStress1", principalStress1)
    call MPAS_pool_get_array(hexdyn, "principalStress2", principalStress2)
    call MPAS_pool_get_array(hexdyn, "replacementPressure", replacementPressure)

    call MPAS_pool_get_array(hexfor, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(hexfor, "vOceanVelocity", vOceanVelocity)

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

    ! calculate the ocean stress
    call ocean_stress(& 
         oceanStressU,   &
         oceanStressV,   &
         uOceanVelocity, & 
         vOceanVelocity, &
         fVertex,        &
         solveVelocity)

    ! calculate the surface tilt force
    call surface_tilt(&
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         uOceanVelocity,    & 
         vOceanVelocity,    &
         totalMassVertex,   &
         fVertex,           &
         solveVelocity)

    ! subcycle the dynamics
    call subcycle_dynamics_hex(mesh, &
                               icestate,   boundary,  &
                               normal,     hexdyn,    &
                               hexweak,    hexvar,    &
                               hexwach,    hexpwl,    &
                               hexfor,                &
                               dtDynamics, dtElastic, &
                               divergenceFormulation, &
                               evpDamping, nstep)

    elapsedTime = elapsedTime + dt

    call principal_stresses_driver(&
         hexweak, &
         hexvar, &
         replacementPressure, &
         divergenceFormulation, &
         principalStress1, &
         principalStress2)

  end subroutine run_dynamics_hex

  !-------------------------------------------------------------

  subroutine principal_stresses_driver(&
       hexweak, &
       hexvar, &
       replacementPressure, &
       divergenceFormulation, &
       principalStress1, &
       principalStress2)

    use cice_dynamics_shared, only: &
         principal_stresses         

    type (MPAS_pool_type), pointer, intent(in) :: &
         hexweak, &
         hexvar

    real(kind=RKIND), dimension(:), intent(in) :: &
         replacementPressure

    character(len=*), intent(in) :: &
         divergenceFormulation
    
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

    ! calculate the principal stresses
    if (trim(divergenceFormulation) == "weak") then

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

    else if (trim(divergenceFormulation) == "variational") then

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
                                  normal,     hexdyn,    &
                                  hexweak,    hexvar,    &
                                  hexwach,    hexpwl,    &
                                  hexfor,                &
                                  dtDynamics, dtElastic, &
                                  divergenceFormulation, &
                                  evpDamping, nstep)

    use cice_dynamics_shared, only: nElasticSubcycle

    type(MPAS_pool_type), pointer, intent(in) :: mesh
    type(MPAS_pool_type), pointer :: icestate
    type(MPAS_pool_type), pointer :: boundary
    type(MPAS_pool_type), pointer :: normal
    type(MPAS_pool_type), pointer :: hexdyn
    type(MPAS_pool_type), pointer :: hexweak
    type(MPAS_pool_type), pointer :: hexvar
    type(MPAS_pool_type), pointer :: hexwach
    type(MPAS_pool_type), pointer :: hexpwl
    type(MPAS_pool_type), pointer :: hexfor

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

    do iElasticSubcycle = 1, nElasticSubcycle

       !write(*,*) "subcycle: ", iElasticSubcycle
       call single_subcycle_dynamics_hex(mesh, &
                                         icestate,   boundary,  &
                                         normal,     hexdyn,    &
                                         hexweak,    hexvar,    &
                                         hexwach,    hexpwl,    &
                                         hexfor,                &
                                         dtDynamics, dtElastic, &
                                         divergenceFormulation, &
                                         evpDamping,            &
                                         iElasticSubcycle, nstep)

    enddo

  end subroutine subcycle_dynamics_hex

  !-------------------------------------------------------------

  subroutine single_subcycle_dynamics_hex(mesh, &
                                          icestate,   boundary,  &
                                          normal,     hexdyn,    &
                                          hexweak,    hexvar,    &
                                          hexwach,    hexpwl,    &
                                          hexfor,                &
                                          dtDynamics, dtElastic, &
                                          divergenceFormulation, &
                                          evpDamping,            &
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
    type(MPAS_pool_type), pointer :: icestate
    type(MPAS_pool_type), pointer :: boundary
    type(MPAS_pool_type), pointer :: normal
    type(MPAS_pool_type), pointer :: hexdyn
    type(MPAS_pool_type), pointer :: hexweak
    type(MPAS_pool_type), pointer :: hexvar
    type(MPAS_pool_type), pointer :: hexwach
    type(MPAS_pool_type), pointer :: hexpwl
    type(MPAS_pool_type), pointer :: hexfor

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

    integer, pointer :: &
         nVertices

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveStress

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaVertex, &
         totalMassVertex, &
         fVertex, &
         latCellRotated, &
         latVertexRotated, &
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
         oceanStressU, &
         oceanStressV, &
         oceanStressCoeff, &
         uOceanVelocity, &
         vOceanVelocity

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         normalVectorPolygon, &
         normalVectorTriangle

    call MPAS_pool_get_array(icestate, "iceAreaVertex", iceAreaVertex)
    call MPAS_pool_get_array(icestate, "totalMassVertex", totalMassVertex)

    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "fvertex", fvertex)

    call MPAS_pool_get_array(normal, "normalVectorPolygon", normalVectorPolygon)
    call MPAS_pool_get_array(normal, "normalVectorTriangle", normalVectorTriangle)
    call MPAS_pool_get_array(normal, "latCellRotated", latCellRotated)
    call MPAS_pool_get_array(normal, "latVertexRotated", latVertexRotated)

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
    call MPAS_pool_get_array(hexdyn, "oceanStressU", oceanStressU)
    call MPAS_pool_get_array(hexdyn, "oceanStressV", oceanStressV)
    call MPAS_pool_get_array(hexdyn, "oceanStressCoeff", oceanStressCoeff)

    call MPAS_pool_get_array(hexfor, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(hexfor, "vOceanVelocity", vOceanVelocity)

    if (trim(divergenceFormulation) == "weak") then

       call internal_stress_weak(&
            mesh,                hexweak, &
            uVelocity,           vVelocity, &
            normalVectorPolygon, normalVectorTriangle, &
            latCellRotated,      latVertexRotated, &
            icePressure,         replacementPressure, &
            dtElastic,           evpDamping, &
            stressDivergenceU,   stressDivergenceV, &
            solveStress,         solveVelocity)

    else if (trim(divergenceFormulation) == "wachspress") then

       call internal_stress_wachspress(&
            mesh, &
            hexvar,              hexwach, &
            uVelocity,           vVelocity, &
            icePressure,         replacementPressure, &
            dtElastic,           evpDamping, &
            stressDivergenceU,   stressDivergenceV, &
            solveStress,         solveVelocity)

    else if (trim(divergenceFormulation) == "pwl") then

       call internal_stress_pwl(&
            mesh, &
            hexvar,              hexpwl, &
            uVelocity,           vVelocity, &
            icePressure,         replacementPressure, &
            dtElastic,           evpDamping, &
            stressDivergenceU,   stressDivergenceV, &
            solveStress,         solveVelocity)

    endif

    call ocean_stress_coefficient(& 
         oceanStressCoeff, &
         uOceanVelocity,   & 
         vOceanVelocity,   &
         uVelocity,        &   
         vVelocity,        &
         iceAreaVertex,    &
         solveVelocity)

    call solve_velocity(&
         nVertices,         &
         solveVelocity,     &
         uVelocity,         &   
         vVelocity,         &
         totalMassVertex,   & 
         fVertex,           &
         stressDivergenceU, & 
         stressDivergenceV, &
         airStressVertexU,  & 
         airStressVertexV,  &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         oceanStressU,      &
         oceanStressV,      &
         oceanStressCoeff,  &
         dtElastic, mesh)

  end subroutine single_subcycle_dynamics_hex

  !-------------------------------------------------------------

  subroutine aggregate_mass_and_area(&
       mesh, &
       iceAreaCategory, iceVolumeCategory, snowVolumeCategory, &
       iceAreaCell, iceVolumeCell, snowVolumeCell, totalMassCell)

    use cice_constants, only: &
         rhoi, &
         rhos

    type (MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(in) :: &
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

    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nCategories", nCategories)

    do iCell = 1, nCells

       totalMassCell(iCell)  = 0.0_RKIND
       iceAreaCell(iCell)    = 0.0_RKIND
       iceVolumeCell(iCell)  = 0.0_RKIND
       snowVolumeCell(iCell) = 0.0_RKIND

       do iCategory = 1, nCategories

          iceAreaCell(iCell)    = iceAreaCell(iCell) + iceAreaCategory(iCategory,iCell)
          iceVolumeCell(iCell)  = iceVolumeCell(iCell) + iceVolumeCategory(iCategory,iCell)
          snowVolumeCell(iCell) = snowVolumeCell(iCell) + snowVolumeCategory(iCategory,iCell)

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

    call MPAS_pool_get_array(mesh, "nCells", nCells)
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

    call MPAS_pool_get_array(mesh, "nVertices", nVertices)

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

