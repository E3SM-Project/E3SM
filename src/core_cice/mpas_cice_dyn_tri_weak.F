module cice_dyn_tri_weak

  use mpas_grid_types

  implicit none

  private
  public :: init_dynamics_tri_weak, &
            strain_tensor_tri_weak, &
            stress_tensor_tri_weak, &
            stress_divergence_tri_weak

contains

  !-------------------------------------------------------------
  ! Initialization
  !-------------------------------------------------------------

  subroutine init_dynamics_tri_weak(mesh, boundary, normal)

    use cice_dynamics_shared, only: interior_vertices, &
                                    normal_vectors

    type(mesh_type), intent(inout) :: mesh

    type(boundary_type), pointer :: boundary
    type(normal_type),   pointer :: normal

    call boundary_cells(mesh, &
         boundary % boundaryCell, &
         boundary % boundaryCell2, &
         boundary % interiorCell) 

    call boundary_vertices(mesh, &
         boundary % boundaryVertex, &
         boundary % interiorVertex, &
         boundary % interiorVertex2, &
         boundary % boundaryCell)

    call normal_vectors(mesh, &
         normal % normalVectorPolygon % array, &
         normal % normalVectorTriangle % array, &
         boundary % interiorVertex % array, &
         "planar")

  end subroutine init_dynamics_tri_weak

  !-------------------------------------------------------------
  ! Time step
  !-------------------------------------------------------------

  subroutine strain_tensor_tri_weak(mesh, &
                                    strain11, strain22, strain12, &
                                    uVelocity, vVelocity, &
                                    normalVectorTriangle, interiorVertex2)

    type(mesh_type), intent(in) :: mesh
    
    real(kind=RKIND), dimension(:), intent(out) :: &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorTriangle

    integer, dimension(:), intent(in) :: &
         interiorVertex2

    integer :: &
         iVertex, &
         iVertexDegree, &
         iEdge, &
         iCellOnEdge, &
         iCell

    real(kind=RKIND) :: &
         uVelocityEdge, &
         vVelocityEdge

    do iVertex = 1, mesh % nVertices
       
       if (interiorVertex2(iVertex) == 1) then
          
          strain11(iVertex) = 0.0_RKIND
          strain22(iVertex) = 0.0_RKIND
          strain12(iVertex) = 0.0_RKIND

          do iVertexDegree = 1, mesh % vertexDegree

             iEdge = mesh % edgesOnVertex % array(iVertexDegree,iVertex)

             uVelocityEdge = 0.0_RKIND
             vVelocityEdge = 0.0_RKIND

             do iCellOnEdge = 1, 2

                iCell = mesh % cellsOnEdge % array(iCellOnEdge,iEdge)

                uVelocityEdge = uVelocityEdge + uVelocity(iCell)
                vVelocityEdge = vVelocityEdge + vVelocity(iCell)

             enddo ! iCellOnEdge

             uVelocityEdge = uVelocityEdge / 2.0_RKIND
             vVelocityEdge = vVelocityEdge / 2.0_RKIND

             strain11(iVertex) = strain11(iVertex) + uVelocityEdge * normalVectorTriangle(1,iVertexDegree,iVertex) * mesh % dcEdge % array(iEdge)
             strain22(iVertex) = strain22(iVertex) + vVelocityEdge * normalVectorTriangle(2,iVertexDegree,iVertex) * mesh % dcEdge % array(iEdge)
             strain12(iVertex) = strain12(iVertex) + 0.5_RKIND * ( &
                  uVelocityEdge * normalVectorTriangle(2,iVertexDegree,iVertex) + &
                  vVelocityEdge * normalVectorTriangle(1,iVertexDegree,iVertex) ) * mesh % dcEdge % array(iEdge)

          enddo ! iVertexDegree

          strain11(iVertex) = strain11(iVertex) / mesh % areaTriangle % array(iVertex)
          strain22(iVertex) = strain22(iVertex) / mesh % areaTriangle % array(iVertex)
          strain12(iVertex) = strain12(iVertex) / mesh % areaTriangle % array(iVertex)
       
       endif ! interiorVertex2

    enddo ! iVertex

  end subroutine strain_tensor_tri_weak

  !-------------------------------------------------------------

  subroutine stress_tensor_tri_weak(mesh, &
                                    stress11, stress22, stress12, &
                                    strain11, strain22, strain12, &
                                    icePressureVertex, &
                                    dtElastic, evpDamping, &
                                    interiorVertex)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         stress11, &
         stress22, &
         stress12

   real(kind=RKIND), dimension(:), intent(in) :: &
         strain11, &
         strain22, &
         strain12, &
         icePressureVertex

    real(kind=RKIND), intent(in) :: &
         dtElastic

    logical, intent(in) :: &
         evpDamping

    integer, dimension(:), intent(in) :: &
         interiorVertex

    integer :: &
         iVertex

    do iVertex = 1, mesh % nVertices

       if (interiorVertex(iVertex) == 1) then

          call evp_constitutive_relation(stress11(iVertex),    stress22(iVertex),    stress12(iVertex), &
                                         strain11(iVertex),    strain22(iVertex),    strain12(iVertex), &
                                         icePressureVertex(iVertex),           &
                                         mesh % areaTriangle % array(iVertex), &
                                         dtElastic, evpDamping)
          
          !call linear_constitutive_relation(stress11(iVertex), stress22(iVertex), stress12(iVertex), &
          !                                  strain11(iVertex), strain22(iVertex), strain12(iVertex))
          
       endif ! interiorVertex
          
    enddo ! iVertex

  end subroutine stress_tensor_tri_weak

  !-------------------------------------------------------------

  subroutine stress_divergence_tri_weak(mesh, &
                                        stressDivergenceU,    stressDivergenceV, &
                                        stress11,             stress22,          &
                                        stress12,                                &
                                        normalVectorPolygon)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, &
         stressDivergenceV

    real(kind=RKIND), dimension(:), intent(in) :: &
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorPolygon

    real(kind=RKIND) :: &
         stress11Edge, &
         stress22Edge, &
         stress12Edge

    integer :: &
         iCell, &
         iEdgeOnCell, &
         iEdge, &
         iVertexOnEdge, &
         iVertex

    do iCell = 1, mesh % nCells

       stressDivergenceU(iCell) = 0.0_RKIND
       stressDivergenceV(iCell) = 0.0_RKIND

       do iEdgeOnCell = 1, mesh % nEdgesOnCell % array(iCell)

          iEdge = mesh % edgesOnCell % array(iEdgeOnCell,iCell)

          stress11Edge = 0.0_RKIND
          stress22Edge = 0.0_RKIND
          stress12Edge = 0.0_RKIND

          do iVertexOnEdge = 1, 2

             iVertex = mesh % verticesOnEdge % array(iVertexOnEdge,iEdge)

             stress11Edge = stress11Edge + stress11(iVertex)
             stress22Edge = stress22Edge + stress22(iVertex)
             stress12Edge = stress12Edge + stress12(iVertex)

          enddo ! iVertexOnEdge

          stress11Edge = stress11Edge / 2.0_RKIND
          stress22Edge = stress22Edge / 2.0_RKIND
          stress12Edge = stress12Edge / 2.0_RKIND

          stressDivergenceU(iCell) = stressDivergenceU(iCell) + &
               (stress11Edge * normalVectorPolygon(1,iEdgeOnCell,iCell) + &
                stress12Edge * normalVectorPolygon(2,iEdgeOnCell,iCell)) * mesh % dvEdge % array(iEdge)

          stressDivergenceV(iCell) = stressDivergenceV(iCell) + &
               (stress22Edge * normalVectorPolygon(2,iEdgeOnCell,iCell) + &
                stress12Edge * normalVectorPolygon(1,iEdgeOnCell,iCell)) * mesh % dvEdge % array(iEdge)

       enddo ! iEdgeOnCell

       stressDivergenceU(iCell) = stressDivergenceU(iCell) / mesh % areaCell % array(iCell)
       stressDivergenceV(iCell) = stressDivergenceV(iCell) / mesh % areaCell % array(iCell)

    enddo ! iCell

  end subroutine stress_divergence_tri_weak

  !-------------------------------------------------------------

end module cice_dyn_tri_weak

