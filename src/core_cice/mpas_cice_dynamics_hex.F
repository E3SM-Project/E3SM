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
                               hexvar,   hexpwl,      &
                               divergenceFormulation, &
                               gridCurvature)

    use cice_dynamics_shared, only: init_evp
    use cice_dyn_hex_weak, only: init_dynamics_hex_weak
    use cice_dyn_hex_variational, only: init_dynamics_hex_variational
    use cice_dyn_hex_pwl, only: init_dynamics_hex_pwl

    use cice_diagnostics, only: plot_pop_grid_vertices

    type(MPAS_pool_type), intent(inout) :: mesh

    type(MPAS_pool_type), pointer :: normal
    type(MPAS_pool_type), pointer :: hexvar
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
                              evpDamping, &
                              nstep)

    use cice_mesh, only: interpolate_cell_to_vertex

    use cice_dynamics_shared, only: ice_strength, &           
                                    air_stress, &
                                    ocean_stress, &
                                    surface_tilt, &
                                    principal_stresses, &
                                    dtDynamics, &
                                    dtElastic, &
                                    elapsedTime

    use cice_diagnostics, only: gnuplot_triangle, gnuplot_cell, writeout_minmax, plot_pop_grid_cells, &
         plot_cells_spherical, &
         plot_vertices_spherical, &
         plot_pop_grid_vertices

    type(MPAS_pool_type), intent(in) :: mesh

    type(MPAS_pool_type), pointer :: icestate
    type(MPAS_pool_type), pointer :: boundary
    type(MPAS_pool_type), pointer :: normal
    type(MPAS_pool_type), pointer :: atmos_forcing    
    type(MPAS_pool_type), pointer :: hexdyn
    type(MPAS_pool_type), pointer :: hexweak
    type(MPAS_pool_type), pointer :: hexvar
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

    real(kind=RKIND), dimension(:), allocatable :: &
         iceAreaCell, &
         iceVolumeCell

    integer, pointer :: &
         nCells, &
         index_iceAreaCell, &
         index_iceVolumeCell

    integer, dimension(:), pointer :: &
         solveStress, &
         solveVelocity, &
         interiorVertex

    real(kind=RKIND), dimension(:), pointer :: &
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
         fVertex, &
         stress11Weak, &
         stress22Weak, &
         stress12Weak

    real(kind=RKIND), dimension(:,:), pointer :: &
         stress11Var, &
         stress22Var, &
         stress12Var

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         tracers

    ! init variables
    call MPAS_pool_get_array(icestate, "tracers", tracers)
    call MPAS_pool_get_array(icestate, "index_iceAreaCell", index_iceAreaCell)
    call MPAS_pool_get_array(icestate, "index_iceVolumeCell", index_iceVolumeCell)
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

    call MPAS_pool_get_array(hexweak, "stress11", stress11Weak) 
    call MPAS_pool_get_array(hexweak, "stress22", stress22Weak) 
    call MPAS_pool_get_array(hexweak, "stress12", stress12Weak) 

    call MPAS_pool_get_array(hexvar, "stress11", stress11Var) 
    call MPAS_pool_get_array(hexvar, "stress22", stress22Var)
    call MPAS_pool_get_array(hexvar, "stress12", stress12Var)

    allocate(iceAreaCell(1:nCells))
    allocate(iceVolumeCell(1:nCells))

    iceAreaCell   = tracers(index_iceAreaCell, 1, :)
    iceVolumeCell = tracers(index_iceVolumeCell, 1, :)

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
                               hexpwl,     hexfor,    &
                               dtDynamics, dtElastic, &
                               divergenceFormulation, &
                               evpDamping, nstep)

    elapsedTime = elapsedTime + dt

    ! calculate the principal stresses
    if (trim(divergenceFormulation) == "weak") then

       call principal_stresses(&
            principalStress1, &
            principalStress2, &
            stress11Weak,     &
            stress22Weak,     &
            stress12Weak,     &
            replacementPressure)

    else if (trim(divergenceFormulation) == "variational") then

       call principal_stresses(&
            principalStress1, &
            principalStress2, &
            stress11Var(1,:), &
            stress22Var(1,:), &
            stress12Var(1,:), &
            replacementPressure)

    endif

    deallocate(iceAreaCell)
    deallocate(iceVolumeCell)

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

    type(MPAS_pool_type), intent(in) :: mesh

    type(MPAS_pool_type), pointer :: icestate
    type(MPAS_pool_type), pointer :: boundary
    type(MPAS_pool_type), pointer :: normal
    type(MPAS_pool_type), pointer :: hexdyn
    type(MPAS_pool_type), pointer :: hexweak
    type(MPAS_pool_type), pointer :: hexvar
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
                                         hexpwl,     hexfor,    &
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
                                          hexpwl,     hexfor,    &
                                          dtDynamics, dtElastic, &
                                          divergenceFormulation, &
                                          evpDamping,            &
                                          iElasticSubcycle, nstep)

    use cice_mesh, only: interpolate_cell_to_vertex                                    

    use cice_dynamics_shared, only: ocean_stress_coefficient, &
                                    solve_velocity

    use cice_dyn_hex_weak, only: strain_tensor_hex_weak, &
                                 stress_tensor_hex_weak, &
                                 stress_divergence_hex_weak

    use cice_dyn_hex_variational, only: strain_tensor_hex_variational, &
                                        stress_tensor_hex_variational, &
                                        stress_divergence_hex_variational

    use cice_dyn_hex_pwl, only: strain_tensor_hex_pwl, &
                                stress_divergence_hex_pwl

    use cice_diagnostics, only: &
         plot_cells_spherical, &
         plot_vertices_spherical, &
         plot_pop_grid_vertices, &
         plot_pop_grid_cells

    type(MPAS_pool_type), intent(in) :: mesh

    type(MPAS_pool_type), pointer :: icestate
    type(MPAS_pool_type), pointer :: boundary
    type(MPAS_pool_type), pointer :: normal
    type(MPAS_pool_type), pointer :: hexdyn
    type(MPAS_pool_type), pointer :: hexweak
    type(MPAS_pool_type), pointer :: hexvar
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

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex

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
         vOceanVelocity, &
         strain11Weak, &
         strain22Weak, &
         strain12Weak, &
         stress11Weak, &
         stress22Weak, &
         stress12Weak

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11Var, &
         strain22Var, &
         strain12Var, &
         stress11Var, &
         stress22Var, &
         stress12Var, &
         basisSubArea

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         normalVectorPolygon, &
         normalVectorTriangle, &
         basisGradientU, &
         basisGradientV, &
         basisIntegralsU, &
         basisIntegralsV, &
         basisGradientUPWL, &
         basisGradientVPWL

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

    call MPAS_pool_get_array(hexweak, "strain11", strain11Weak)
    call MPAS_pool_get_array(hexweak, "strain22", strain22Weak)
    call MPAS_pool_get_array(hexweak, "strain12", strain12Weak)
    call MPAS_pool_get_array(hexweak, "stress11", stress11Weak)
    call MPAS_pool_get_array(hexweak, "stress22", stress22Weak)
    call MPAS_pool_get_array(hexweak, "stress12", stress12Weak)

    call MPAS_pool_get_array(hexvar, "strain11", strain11Var)
    call MPAS_pool_get_array(hexvar, "strain22", strain22Var)
    call MPAS_pool_get_array(hexvar, "strain12", strain12Var)
    call MPAS_pool_get_array(hexvar, "stress11", stress11Var)
    call MPAS_pool_get_array(hexvar, "stress22", stress22Var)
    call MPAS_pool_get_array(hexvar, "stress12", stress12Var)
    call MPAS_pool_get_array(hexvar, "basisGradientU", basisGradientU)
    call MPAS_pool_get_array(hexvar, "basisGradientV", basisGradientV)
    call MPAS_pool_get_array(hexvar, "basisIntegralsU", basisIntegralsU)
    call MPAS_pool_get_array(hexvar, "basisIntegralsV", basisIntegralsV)
    call MPAS_pool_get_array(hexvar, "cellVerticesAtVertex", cellVerticesAtVertex)

    call MPAS_pool_get_array(hexpwl, "basisGradientU", basisGradientUPWL)
    call MPAS_pool_get_array(hexpwl, "basisGradientV", basisGradientVPWL)
    call MPAS_pool_get_array(hexpwl, "basisSubArea", basisSubArea)

    if (trim(divergenceFormulation) == "weak") then

       !call strain_tensor_hex_weak(mesh, &
       !     strain11Weak,        &
       !     strain22Weak,        &
       !     strain12Weak,        &
       !     uVelocity,           &
       !     vVelocity,           &
       !     normalVectorPolygon, &
       !     latCellRotated,      &
       !     solveStress)

       !call stress_tensor_hex_weak(mesh, &
       !     stress11Weak,         &
       !     stress22Weak,         &
       !     stress12Weak,         &
       !     strain11Weak,         &
       !     strain22Weak,         &
       !     strain12Weak,         &
       !     icePressure,          &
       !     replacementPressure,  &
       !     solveStress,          &       
       !     dtElastic, evpDamping)

       !call stress_divergence_hex_weak(mesh, &
       !     stressDivergenceU,    &
       !     stressDivergenceV,    &
       !     stress11Weak,         &
       !     stress22Weak,         &
       !     stress12Weak,         &
       !     normalVectorTriangle, &
       !     latVertexRotated,     &
       !     solveVelocity)

    else if (trim(divergenceFormulation) == "variational") then

       call strain_tensor_hex_variational(mesh, &
            strain11Var,    &
            strain22Var,    &
            strain12Var,    &
            uVelocity,      &
            vVelocity,      &
            basisGradientU, &
            basisGradientV, &
            solveStress)

       call stress_tensor_hex_variational(mesh, &
            stress11Var,          &
            stress22Var,          &
            stress12Var,          &
            strain11Var,          &
            strain22Var,          &
            strain12Var,          &
            icePressure,          &
            replacementPressure,  &
            solveStress,          &
            dtElastic, evpDamping)

       call stress_divergence_hex_variational(mesh, &
            stressDivergenceU,    &
            stressDivergenceV,    &
            stress11Var,          &
            stress22Var,          &
            stress12Var,          &
            basisIntegralsU,      &
            basisIntegralsV,      &
            cellVerticesAtVertex, &
            solveVelocity)

    else if (trim(divergenceFormulation) == "pwl") then

       call strain_tensor_hex_pwl(mesh, &
            strain11Var,       &
            strain22Var,       &
            strain12Var,       &
            uVelocity,         &
            vVelocity,         &
            basisGradientUPWL, &
            basisGradientVPWL, &
            solveStress)

       call stress_tensor_hex_variational(mesh, &
            stress11Var,          &
            stress22Var,          &
            stress12Var,          &
            strain11Var,          &
            strain22Var,          &
            strain12Var,          &
            icePressure,          &
            replacementPressure,  &
            solveStress,          &
            dtElastic, evpDamping)

       call stress_divergence_hex_pwl(mesh, &
            stressDivergenceU,    &
            stressDivergenceV,    &
            stress11Var,          &
            stress22Var,          &
            stress12Var,          &
            basisGradientUPWL,    &
            basisGradientVPWL,    &
            basisSubArea,         &
            cellVerticesAtVertex, &
            solveVelocity)

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

  subroutine stress_calculation_mask_hex(mesh,     solveStress, &
                                         areaCell, massCell)

    use cice_dynamics_shared, only: areaMinimum, massMinimum

    type(MPAS_pool_type), intent(in) :: mesh

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

    type(MPAS_pool_type), intent(in) :: mesh

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

