module cice_dyn_hex_weak

  use mpas_grid_types

  implicit none

  private
  public :: &
       init_dynamics_hex_weak, &
       internal_stress_weak, &
       strain_tensor_hex_weak, &
       stress_tensor_hex_weak, &
       stress_divergence_hex_weak

contains

  !-------------------------------------------------------------
  ! Initialization
  !-------------------------------------------------------------

  subroutine init_dynamics_hex_weak(mesh, boundary, normal, gridCurvature)

    use cice_mesh, only: normal_vectors

    type(MPAS_pool_type), pointer, intent(inout) :: mesh

    type(MPAS_pool_type), pointer :: normal
    type(MPAS_pool_type), pointer :: boundary

    character(len=*), intent(in) :: gridCurvature

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         normalVectorPolygon, &
         normalVectorTriangle

    real(kind=RKIND), dimension(:), pointer :: &
         latCellRotated, &
         latVertexRotated

    integer, dimension(:), pointer :: &
         interiorVertex

    call MPAS_pool_get_array(mesh, "normalVectorPolygon", normalVectorPolygon)
    call MPAS_pool_get_array(mesh, "normalVectorTriangle", normalVectorTriangle)
    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call normal_vectors(mesh, &
         normalVectorPolygon,         &
         normalVectorTriangle,        &
         latCellRotated,              &
         latVertexRotated,            &
         interiorVertex,            &
         gridCurvature)

  end subroutine init_dynamics_hex_weak

  !-------------------------------------------------------------
  ! Time step
  !-------------------------------------------------------------

  subroutine internal_stress_weak(&
       mesh,                hexweak, &
       uVelocity,           vVelocity, &
       normalVectorPolygon, normalVectorTriangle, &
       latCellRotated,      latVertexRotated, &
       icePressure,         replacementPressure, &
       dtElastic,           evpDamping, &
       stressDivergenceU,   stressDivergenceV, &
       solveStress,         solveVelocity)

    type (MPAS_pool_type), pointer :: &
         mesh, &
         hexweak

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, &
         vVelocity, &
         latCellRotated, &
         latVertexRotated         

    real(kind=RKIND), dimension(:), intent(inout) :: &
         icePressure, &
         replacementPressure

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorPolygon, &
         normalVectorTriangle

    real(kind=RKIND), intent(in) :: &
         dtElastic

    logical, intent(in) :: &
         evpDamping

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, &
         stressDivergenceV

    integer, dimension(:), intent(in) :: &
         solveStress, &
         solveVelocity

    real(kind=RKIND), dimension(:), pointer :: &
         strain11, &
         strain22, &
         strain12, &
         stress11, &
         stress22, &
         stress12

    call MPAS_pool_get_array(hexweak, "strain11", strain11)
    call MPAS_pool_get_array(hexweak, "strain22", strain22)
    call MPAS_pool_get_array(hexweak, "strain12", strain12)
    call MPAS_pool_get_array(hexweak, "stress11", stress11)
    call MPAS_pool_get_array(hexweak, "stress22", stress22)
    call MPAS_pool_get_array(hexweak, "stress12", stress12)

    call strain_tensor_hex_weak(mesh, &
         strain11,            &
         strain22,            &
         strain12,            &
         uVelocity,           &
         vVelocity,           &
         normalVectorPolygon, &
         latCellRotated,      &
         solveStress)
    
    call stress_tensor_hex_weak(mesh, &
         stress11,             &
         stress22,             &
         stress12,             &
         strain11,             &
         strain22,             &
         strain12,             &
         icePressure,          &
         replacementPressure,  &
         solveStress,          &       
         dtElastic, evpDamping)

    call stress_divergence_hex_weak(mesh, &
         stressDivergenceU,    &
         stressDivergenceV,    &
         stress11,             &
         stress22,             &
         stress12,             &
         normalVectorTriangle, &
         latVertexRotated,     &
         solveVelocity)

  end subroutine internal_stress_weak

  !-------------------------------------------------------------

  subroutine strain_tensor_hex_weak(mesh, &
                                    strain11, strain22,   &
                                    strain12,             &
                                    uVelocity, vVelocity, &
                                    normalVectorPolygon,  &
                                    latCellRotated,       &
                                    solveStress)

    use cice_constants, only: deg_to_rad

    type(MPAS_pool_type), pointer, intent(in) :: mesh

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

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell, &
         POPindxi, &
         POPindxj

    integer, dimension(:,:), pointer :: &
         verticesOnCell, &
         edgesOnCell, &
         verticesOnEdge

    real(kind=RKIND), pointer :: &
         sphere_radius

    real(kind=RKIND), dimension(:), pointer :: &
         dvEdge, &
         areaCell

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "POPindxi", POPindxi)
    call MPAS_pool_get_array(mesh, "POPindxj", POPindxj)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnEdge", verticesOnEdge)
    call MPAS_pool_get_array(mesh, "sphere_radius", sphere_radius)
    call MPAS_pool_get_array(mesh, "dvEdge", dvEdge)
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    do iCell = 1, nCells

       strain11(iCell) = 0.0_RKIND
       strain22(iCell) = 0.0_RKIND
       strain12(iCell) = 0.0_RKIND

       if (solveStress(iCell) == 1) then

          uCellCentre = 0.0_RKIND
          vCellCentre = 0.0_RKIND

          do iEdgeOnCell = 1, nEdgesOnCell(iCell)

             ! cell centre velocities
             iVertex = verticesOnCell(iEdgeOnCell,iCell)

             uCellCentre = uCellCentre + uVelocity(iVertex)
             vCellCentre = vCellCentre + vVelocity(iVertex)

             ! interpolated edge velocity
             iEdge = edgesOnCell(iEdgeOnCell,iCell)

             uVelocityEdge = 0.0_RKIND
             vVelocityEdge = 0.0_RKIND

             do iVertexOnEdge = 1, 2

                iVertex = verticesOnEdge(iVertexOnEdge,iEdge)

                uVelocityEdge = uVelocityEdge + uVelocity(iVertex)
                vVelocityEdge = vVelocityEdge + vVelocity(iVertex)

             enddo ! iVertexOnEdge

             uVelocityEdge = uVelocityEdge / 2.0_RKIND
             vVelocityEdge = vVelocityEdge / 2.0_RKIND

             ! summation over edges
             strain11(iCell) = strain11(iCell) + uVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) * dvEdge(iEdge)
             strain22(iCell) = strain22(iCell) + vVelocityEdge * normalVectorPolygon(2,iEdgeOnCell,iCell) * dvEdge(iEdge)
             strain12(iCell) = strain12(iCell) + 0.5_RKIND * ( &
                  uVelocityEdge * normalVectorPolygon(2,iEdgeOnCell,iCell) + &
                  vVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) ) * dvEdge(iEdge)

             !write(72,*) POPindxi(iCell), POPindxj(iCell), iEdgeOnCell, &
             !     uVelocityEdge, normalVectorPolygon(1,iEdgeOnCell,iCell), dvEdge(iEdge), &
             !     uVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) * dvEdge(iEdge)

          enddo ! iEdgeOnCell

          uCellCentre = uCellCentre / real(nEdgesOnCell(iCell), RKIND)
          vCellCentre = vCellCentre / real(nEdgesOnCell(iCell), RKIND)

          !write(72,*) POPindxi(iCell), POPindxj(iCell), 0, &
          !     strain11(iCell), strain11(iCell) / areaCell(iCell), areaCell(iCell)

          strain11(iCell) = strain11(iCell) / areaCell(iCell)
          strain22(iCell) = strain22(iCell) / areaCell(iCell)
          strain12(iCell) = strain12(iCell) / areaCell(iCell)

          ! metric terms
          strain11(iCell) = strain11(iCell) - (vCellCentre * tan(latCellRotated(iCell))) / sphere_radius
          strain12(iCell) = strain12(iCell) + (uCellCentre * tan(latCellRotated(iCell)) * 0.5_RKIND) / sphere_radius

          if (abs(strain11(iCell)) < 1e-10_RKIND) strain11(iCell) = 0.0_RKIND
          write(71,*) POPindxi(iCell), POPindxj(iCell), strain11(iCell), strain22(iCell), strain12(iCell)

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

    type(MPAS_pool_type), pointer, intent(in) :: mesh

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

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    do iCell = 1, nCells

       if (solveStress(iCell) == 1) then

          call evp_constitutive_relation(stress11(iCell),    stress22(iCell),    stress12(iCell), &
                                         strain11(iCell),    strain22(iCell),    strain12(iCell), &
                                         icePressure(iCell), replacementPressure(iCell),          &
                                         areaCell(iCell),                                         &
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

    type(MPAS_pool_type), pointer, intent(in) :: mesh

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

    integer, pointer :: &
         nVertices, &
         vertexDegree

    integer, dimension(:), pointer :: &
         POPindxiv, &
         POPindxjv

    integer, dimension(:,:), pointer :: &
         cellsOnVertex, &
         edgesOnVertex, &
         cellsOnEdge

    real(kind=RKIND), pointer :: &
         sphere_radius
    
    real(kind=RKIND), dimension(:), pointer :: &
         areaTriangle, &
         dcEdge

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)

    call MPAS_pool_get_array(mesh, "POPindxiv", POPindxiv)
    call MPAS_pool_get_array(mesh, "POPindxjv", POPindxjv)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "edgesOnVertex", edgesOnVertex)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "sphere_radius", sphere_radius)
    call MPAS_pool_get_array(mesh, "areaTriangle", areaTriangle)
    call MPAS_pool_get_array(mesh, "dcEdge", dcEdge)

    do iVertex = 1, nVertices

       stressDivergenceU(iVertex) = 0.0_RKIND
       stressDivergenceV(iVertex) = 0.0_RKIND

       if (solveVelocity(iVertex) == 1) then

          stress11Vertex = 0.0_RKIND
          stress22Vertex = 0.0_RKIND
          stress12Vertex = 0.0_RKIND
          
          do iVertexDegree = 1, vertexDegree
             
             ! vertex stresses
             iCell = cellsOnVertex(iVertexDegree,iVertex)

             stress11Vertex = stress11Vertex + stress11(iCell)
             stress22Vertex = stress22Vertex + stress22(iCell)
             stress12Vertex = stress12Vertex + stress12(iCell)
             
             ! interpolated edge velocity
             iEdge = edgesOnVertex(iVertexDegree,iVertex)
             
             stress11Edge = 0.0_RKIND
             stress22Edge = 0.0_RKIND
             stress12Edge = 0.0_RKIND
             
             do iCellOnEdge = 1, 2
                
                iCell = cellsOnEdge(iCellOnEdge,iEdge)
                
                stress11Edge = stress11Edge + stress11(iCell)
                stress22Edge = stress22Edge + stress22(iCell)
                stress12Edge = stress12Edge + stress12(iCell)
                
             enddo ! iCellOnEdge
             
             stress11Edge = stress11Edge / 2.0_RKIND
             stress22Edge = stress22Edge / 2.0_RKIND
             stress12Edge = stress12Edge / 2.0_RKIND
             
             stressDivergenceU(iVertex) = stressDivergenceU(iVertex) + &
                  (stress11Edge * normalVectorTriangle(1,iVertexDegree,iVertex) + &
                   stress12Edge * normalVectorTriangle(2,iVertexDegree,iVertex)) * dcEdge(iEdge)
             
             stressDivergenceV(iVertex) = stressDivergenceV(iVertex) + &
                  (stress22Edge * normalVectorTriangle(2,iVertexDegree,iVertex) + &
                   stress12Edge * normalVectorTriangle(1,iVertexDegree,iVertex)) * dcEdge(iEdge)

             write(70,*) POPindxiv(iVertex), POPindxjv(iVertex), &
                  -1, iVertexDegree, stress11Edge, normalVectorTriangle(1,iVertexDegree,iVertex), stress12Edge, normalVectorTriangle(2,iVertexDegree,iVertex), dcEdge(iEdge)

          enddo ! iVertexDegree

          stress11Vertex = stress11Vertex / real(vertexDegree, RKIND)
          stress22Vertex = stress22Vertex / real(vertexDegree, RKIND)
          stress12Vertex = stress12Vertex / real(vertexDegree, RKIND)
          
          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) / areaTriangle(iVertex)
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) / areaTriangle(iVertex)

          write(70,*) POPindxiv(iVertex), POPindxjv(iVertex), &
               0, stressDivergenceU(iVertex)

          ! metric terms
          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) - &
               (tan(latVertexRotated(iVertex)) * stress12Vertex * 2.0_RKIND) / sphere_radius
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) + &
               (tan(latVertexRotated(iVertex)) * (stress11Vertex - stress22Vertex)) / sphere_radius

          write(70,*) POPindxiv(iVertex), POPindxjv(iVertex), &
               1, stressDivergenceU(iVertex)

       endif ! solveVelocity

    enddo ! iVertex

  end subroutine stress_divergence_hex_weak

  !-------------------------------------------------------------

end module cice_dyn_hex_weak


