module cice_dyn_hex_weak

  use mpas_grid_types

  implicit none

  private
  public :: init_dynamics_hex_weak, &
            strain_tensor_hex_weak, &
            stress_tensor_hex_weak, &
            stress_divergence_hex_weak

contains

  !-------------------------------------------------------------
  ! Initialization
  !-------------------------------------------------------------

  subroutine init_dynamics_hex_weak(mesh, boundary, normal, gridCurvature)

    use cice_mesh, only: normal_vectors

    type(mesh_type), intent(inout) :: mesh

    type(normal_type),   pointer :: normal
    type(boundary_type), pointer :: boundary

    character(len=*), intent(in) :: gridCurvature

    call normal_vectors(mesh, &
         normal % normalVectorPolygon % array,         &
         normal % normalVectorTriangle % array,        &
         normal % latCellRotated % array,              &
         normal % latVertexRotated % array,            &
         boundary % interiorVertex % array,            &
         gridCurvature)

  end subroutine init_dynamics_hex_weak

  !-------------------------------------------------------------
  ! Time step
  !-------------------------------------------------------------

  subroutine strain_tensor_hex_weak(mesh, &
                                    strain11, strain22,   &
                                    strain12,             &
                                    uVelocity, vVelocity, &
                                    normalVectorPolygon,  &
                                    latCellRotated,       &
                                    solveStress)

    use cice_constants, only: deg_to_rad

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, &
         vVelocity, &
         latCellRotated

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorPolygon

    integer, dimension(:), intent(in) :: &
         solveStress

    integer :: &
         iCell, &
         iEdgeOnCell, &
         iEdge, &
         iVertexOnEdge, &
         iVertex
    
    real(kind=RKIND) :: &
         uVelocityEdge, &
         vVelocityEdge, &
         uCellCentre, &
         vCellCentre
     
    do iCell = 1, mesh % nCells

       strain11(iCell) = 0.0_RKIND
       strain22(iCell) = 0.0_RKIND
       strain12(iCell) = 0.0_RKIND

       if (solveStress(iCell) == 1) then

          uCellCentre = 0.0_RKIND
          vCellCentre = 0.0_RKIND

          do iEdgeOnCell = 1, mesh % nEdgesOnCell % array(iCell)

             ! cell centre velocities
             iVertex = mesh % verticesOnCell % array(iEdgeOnCell,iCell)

             uCellCentre = uCellCentre + uVelocity(iVertex)
             vCellCentre = vCellCentre + vVelocity(iVertex)

             ! interpolated edge velocity
             iEdge = mesh % edgesOnCell % array(iEdgeOnCell,iCell)

             uVelocityEdge = 0.0_RKIND
             vVelocityEdge = 0.0_RKIND

             do iVertexOnEdge = 1, 2

                iVertex = mesh % verticesOnEdge % array(iVertexOnEdge,iEdge)

                uVelocityEdge = uVelocityEdge + uVelocity(iVertex)
                vVelocityEdge = vVelocityEdge + vVelocity(iVertex)

             enddo ! iVertexOnEdge

             uVelocityEdge = uVelocityEdge / 2.0_RKIND
             vVelocityEdge = vVelocityEdge / 2.0_RKIND

             ! summation over edges
             strain11(iCell) = strain11(iCell) + uVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) * mesh % dvEdge % array(iEdge)
             strain22(iCell) = strain22(iCell) + vVelocityEdge * normalVectorPolygon(2,iEdgeOnCell,iCell) * mesh % dvEdge % array(iEdge)
             strain12(iCell) = strain12(iCell) + 0.5_RKIND * ( &
                  uVelocityEdge * normalVectorPolygon(2,iEdgeOnCell,iCell) + &
                  vVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) ) * mesh % dvEdge % array(iEdge)

             !write(72,*) mesh % POPindxi % array(iCell), mesh % POPindxj % array(iCell), iEdgeOnCell, &
             !     uVelocityEdge, normalVectorPolygon(1,iEdgeOnCell,iCell), mesh % dvEdge % array(iEdge), &
             !     uVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) * mesh % dvEdge % array(iEdge)

          enddo ! iEdgeOnCell

          uCellCentre = uCellCentre / real(mesh % nEdgesOnCell % array(iCell), RKIND)
          vCellCentre = vCellCentre / real(mesh % nEdgesOnCell % array(iCell), RKIND)

          !write(72,*) mesh % POPindxi % array(iCell), mesh % POPindxj % array(iCell), 0, &
          !     strain11(iCell), strain11(iCell) / mesh % areaCell % array(iCell), mesh % areaCell % array(iCell)

          strain11(iCell) = strain11(iCell) / mesh % areaCell % array(iCell)
          strain22(iCell) = strain22(iCell) / mesh % areaCell % array(iCell)
          strain12(iCell) = strain12(iCell) / mesh % areaCell % array(iCell)

          ! metric terms
          strain11(iCell) = strain11(iCell) - (vCellCentre * tan(latCellRotated(iCell))) / mesh % sphere_radius
          strain12(iCell) = strain12(iCell) + (uCellCentre * tan(latCellRotated(iCell)) * 0.5_RKIND) / mesh % sphere_radius

          if (abs(strain11(iCell)) < 1e-10_RKIND) strain11(iCell) = 0.0_RKIND
          write(71,*) mesh % POPindxi % array(iCell), mesh % POPindxj % array(iCell), strain11(iCell), strain22(iCell), strain12(iCell)

       endif ! solveStress

    enddo ! iCell

  end subroutine strain_tensor_hex_weak

  !-------------------------------------------------------------

  subroutine stress_tensor_hex_weak(mesh, &
                                    stress11, stress22, stress12,     &
                                    strain11, strain22, strain12,     &
                                    icePressure, replacementPressure, &
                                    solveStress,                      &
                                    dtElastic,   evpDamping)

    use cice_dynamics_shared, only: evp_constitutive_relation, &
                                    linear_constitutive_relation

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(inout) :: &
         stress11, &
         stress22, &
         stress12, &
         replacementPressure

    real(kind=RKIND), dimension(:), intent(inout) :: &
         strain11, &
         strain22, &
         strain12, &
         icePressure

    integer, dimension(:), intent(in) :: &
         solveStress

    real(kind=RKIND), intent(in) :: &
         dtElastic

    logical, intent(in) :: &
         evpDamping

    integer :: iCell

    do iCell = 1, mesh % nCells

       if (solveStress(iCell) == 1) then

          call evp_constitutive_relation(stress11(iCell),    stress22(iCell),    stress12(iCell), &
                                         strain11(iCell),    strain22(iCell),    strain12(iCell), &
                                         icePressure(iCell), replacementPressure(iCell),          &
                                         mesh % areaCell % array(iCell),                          &
                                         dtElastic,          evpDamping)

          !call linear_constitutive_relation(stress11(iCell), stress22(iCell), stress12(iCell), &
          !                                  strain11(iCell), strain22(iCell), strain12(iCell))

       else

          stress11(iCell) = 0.0_RKIND
          stress22(iCell) = 0.0_RKIND
          stress12(iCell) = 0.0_RKIND

       endif ! solveStress

    end do ! iCell

  end subroutine stress_tensor_hex_weak
  
  !-------------------------------------------------------------

  subroutine stress_divergence_hex_weak(mesh, &
                                        stressDivergenceU,    stressDivergenceV, &
                                        stress11,             stress22,          &
                                        stress12,                                &
                                        normalVectorTriangle, latVertexRotated,  &
                                        solveVelocity)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, &
         stressDivergenceV

    real(kind=RKIND), dimension(:), intent(in) :: &
         stress11, &
         stress22, &
         stress12, &
         latVertexRotated

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorTriangle

    integer, dimension(:), intent(in) :: &
         solveVelocity

    real(kind=RKIND) :: &
         stress11Edge, &
         stress22Edge, &
         stress12Edge, &
         stress11Vertex, &
         stress22Vertex, &
         stress12Vertex

    integer :: &
         iVertex, &
         iVertexDegree, &
         iEdge, &
         iCellOnEdge, &
         iCell

    do iVertex = 1, mesh % nVertices

       stressDivergenceU(iVertex) = 0.0_RKIND
       stressDivergenceV(iVertex) = 0.0_RKIND

       if (solveVelocity(iVertex) == 1) then

          stress11Vertex = 0.0_RKIND
          stress22Vertex = 0.0_RKIND
          stress12Vertex = 0.0_RKIND
          
          do iVertexDegree = 1, mesh % vertexDegree
             
             ! vertex stresses
             iCell = mesh % cellsOnVertex % array(iVertexDegree,iVertex)

             stress11Vertex = stress11Vertex + stress11(iCell)
             stress22Vertex = stress22Vertex + stress22(iCell)
             stress12Vertex = stress12Vertex + stress12(iCell)
             
             ! interpolated edge velocity
             iEdge = mesh % edgesOnVertex % array(iVertexDegree,iVertex)
             
             stress11Edge = 0.0_RKIND
             stress22Edge = 0.0_RKIND
             stress12Edge = 0.0_RKIND
             
             do iCellOnEdge = 1, 2
                
                iCell = mesh % cellsOnEdge % array(iCellOnEdge,iEdge)
                
                stress11Edge = stress11Edge + stress11(iCell)
                stress22Edge = stress22Edge + stress22(iCell)
                stress12Edge = stress12Edge + stress12(iCell)
                
             enddo ! iCellOnEdge
             
             stress11Edge = stress11Edge / 2.0_RKIND
             stress22Edge = stress22Edge / 2.0_RKIND
             stress12Edge = stress12Edge / 2.0_RKIND
             
             stressDivergenceU(iVertex) = stressDivergenceU(iVertex) + &
                  (stress11Edge * normalVectorTriangle(1,iVertexDegree,iVertex) + &
                   stress12Edge * normalVectorTriangle(2,iVertexDegree,iVertex)) * mesh % dcEdge % array(iEdge)
             
             stressDivergenceV(iVertex) = stressDivergenceV(iVertex) + &
                  (stress22Edge * normalVectorTriangle(2,iVertexDegree,iVertex) + &
                   stress12Edge * normalVectorTriangle(1,iVertexDegree,iVertex)) * mesh % dcEdge % array(iEdge)

             write(70,*) mesh % POPindxiv % array(iVertex), mesh % POPindxjv % array(iVertex), &
                  -1, iVertexDegree, stress11Edge, normalVectorTriangle(1,iVertexDegree,iVertex), stress12Edge, normalVectorTriangle(2,iVertexDegree,iVertex), mesh % dcEdge % array(iEdge)

          enddo ! iVertexDegree

          stress11Vertex = stress11Vertex / real(mesh % vertexDegree, RKIND)
          stress22Vertex = stress22Vertex / real(mesh % vertexDegree, RKIND)
          stress12Vertex = stress12Vertex / real(mesh % vertexDegree, RKIND)
          
          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) / mesh % areaTriangle % array(iVertex)
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) / mesh % areaTriangle % array(iVertex)

          write(70,*) mesh % POPindxiv % array(iVertex), mesh % POPindxjv % array(iVertex), &
               0, stressDivergenceU(iVertex)

          ! metric terms
          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) - &
               (tan(latVertexRotated(iVertex)) * stress12Vertex * 2.0_RKIND) / mesh % sphere_radius
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) + &
               (tan(latVertexRotated(iVertex)) * (stress11Vertex - stress22Vertex)) / mesh % sphere_radius

          write(70,*) mesh % POPindxiv % array(iVertex), mesh % POPindxjv % array(iVertex), &
               1, stressDivergenceU(iVertex)

       endif ! solveVelocity

    enddo ! iVertex

  end subroutine stress_divergence_hex_weak

  !-------------------------------------------------------------

end module cice_dyn_hex_weak


