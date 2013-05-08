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

    use cice_dynamics_shared, only: interior_vertices, &
                                    normal_vectors

    type(mesh_type), intent(inout) :: mesh

    type(normal_type),   pointer :: normal
    type(boundary_type), pointer :: boundary

    character(len=*), intent(in) :: gridCurvature

    call interior_vertices(mesh, &
         boundary % interiorVertex % array)

    call normal_vectors(mesh, &
         normal % normalVectorPolygon % array, &
         normal % normalVectorTriangle % array, &
         boundary % interiorVertex % array, &
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
                                    solveStress)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, &
         vVelocity

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
         vVelocityEdge
    
    do iCell = 1, mesh % nCells

       strain11(iCell) = 0.0_RKIND
       strain22(iCell) = 0.0_RKIND
       strain12(iCell) = 0.0_RKIND

       if (solveStress(iCell) == 1) then

          do iEdgeOnCell = 1, mesh % nEdgesOnCell % array(iCell)

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

             strain11(iCell) = strain11(iCell) + uVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) * mesh % dvEdge % array(iEdge)
             strain22(iCell) = strain22(iCell) + vVelocityEdge * normalVectorPolygon(2,iEdgeOnCell,iCell) * mesh % dvEdge % array(iEdge)
             strain12(iCell) = strain12(iCell) + 0.5_RKIND * ( &
                  uVelocityEdge * normalVectorPolygon(2,iEdgeOnCell,iCell) + &
                  vVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) ) * mesh % dvEdge % array(iEdge)

          enddo ! iEdgeOnCell

          strain11(iCell) = strain11(iCell) / mesh % areaCell % array(iCell)
          strain22(iCell) = strain22(iCell) / mesh % areaCell % array(iCell)
          strain12(iCell) = strain12(iCell) / mesh % areaCell % array(iCell)

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

    real(kind=RKIND), dimension(:), intent(in) :: &
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
                                        normalVectorTriangle, solveVelocity)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, &
         stressDivergenceV

    real(kind=RKIND), dimension(:), intent(in) :: &
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorTriangle

    integer, dimension(:), intent(in) :: &
         solveVelocity

    real(kind=RKIND) :: &
         stress11Edge, &
         stress22Edge, &
         stress12Edge

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
          
          do iVertexDegree = 1, mesh % vertexDegree
             
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

          enddo ! iVertexDegree
          
          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) / mesh % areaTriangle % array(iVertex)
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) / mesh % areaTriangle % array(iVertex)

       endif ! solveVelocity

    enddo ! iVertex

  end subroutine stress_divergence_hex_weak

  !-------------------------------------------------------------

end module cice_dyn_hex_weak


