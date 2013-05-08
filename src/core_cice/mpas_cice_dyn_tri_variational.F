module cice_dynamics_evp_tri

  use mpas_grid_types

  implicit none

  private
  public :: init_dynamics_evp, &
            test_tri_divergence_stress_tensor

contains



  !-------------------------------------------------------------
  ! Initialization
  !-------------------------------------------------------------

  subroutine init_dynamics_evp(mesh, basis, boundary)

    use cice_testing, only: gnuplot_cell, gnuplot_triangle, plot_boundary_triangles, &
                            boundary_locations, labels_verticesOnCell, labels_verticesDegreeOnVertex, &
                            writeout_minmax
    
    use cice_dynamics_shared, only: boundary_cells, boundary_vertices, boundary_edges

    type(mesh_type), intent(in) :: mesh

    type(basis_type),    pointer :: basis
    type(boundary_type), pointer :: boundary

    integer :: iCell, iVertex
   
    call boundary_cells(mesh, &
         boundary % boundaryCell % array, &
         boundary % boundaryCell2 % array) 

    call boundary_vertices(mesh, &
         boundary % boundaryVertex % array, &
         boundary % interiorVertex % array, &
         boundary % interiorVertex2 % array, &
         boundary % boundaryCell % array)

    call boundary_edges(mesh, &
         boundary % boundaryEdge % array)

    !call gnuplot_cell(mesh, mesh % xCell % array, "cells.txt", .false., .true.)
    call gnuplot_triangle(mesh, mesh % xVertex % array+mesh % yVertex % array, "triangles.txt", .false., .true.)

    call writeout_minmax()

    !call labels_verticesOnCell(mesh, mesh % verticesOnCell % array, "cell_vertices.txt") 
    !call labels_verticesDegreeOnVertex(mesh, mesh % cellsOnVertex % array, boundary % interiorVertex % array, "triangle_corners.txt") 

    call calc_triangle_corner_at_cell_center(mesh, basis % triangleCornerAtCellCenter % array)
    call calc_edge_opposite_triangle_corner(mesh, basis % edgeOppositeTriangleCorner % array, boundary % interiorVertex % array)

    !call labels_verticesOnCell(mesh, basis % triangleCornerAtCellCenter % array, "triangleCornerAtCellCenter.txt") 
    !call labels_verticesDegreeOnVertex(mesh, basis % edgeOppositeTriangleCorner % array, boundary % interiorVertex % array, "edgeOppositeTriangleCorner.txt") 


    !stop


    call basis_gradients_planar(mesh, &
         basis % basisGradientU % array, &
         basis % basisGradientV % array, &
         boundary % interiorVertex % array)

    call basis_integrals(mesh, &
         basis % basisIntegrals % array, &
         basis % edgeOppositeTriangleCorner % array)

    ! plotting stuff
    !call gnuplot_cell(mesh, mesh % xCell % array, "xCell.txt")
    !call gnuplot_triangle(mesh, mesh % xVertex % array, "xVertex.txt")

    !call boundary_locations(mesh, &
    !     boundary % boundaryCell % array, &
    !     boundary % boundaryVertex % array, &
    !     boundary % interiorVertex % array, &
    !     boundary % boundaryEdge % array)

    !call plot_boundary_triangles(mesh, &
    !     boundary % boundaryVertex % array, &
    !     boundary % interiorVertex % array, &
    !     boundary % boundaryEdge % array)

  end subroutine init_dynamics_evp

  !-------------------------------------------------------------

  subroutine init_dynamics_hex_variational(mesh, trivar, boundary)

    type(mesh_type), intent(in) :: mesh
    type(trivar_type),   pointer :: trivar
    type(boundary_type), pointer :: boundary

    call calc_triangle_corner_at_cell_center(mesh, &
         trivar % triangleCornerAtCellCenter % array)

    call calc_edge_opposite_triangle_corner(mesh, &
         trivar % edgeOppositeTriangleCorner % array, &
         boundary % interiorVertex % array)

    call basis_gradients_planar(mesh, &
         trivar % basisGradientU % array, &
         trivar % basisGradientV % array, &
         boundary % interiorVertex % array)

    call basis_integrals(mesh, &
         trivar % basisIntegral % array, &
         trivar % edgeOppositeTriangleCorner % array)

  end subroutine init_dynamics_hex_variational

  !-------------------------------------------------------------

  subroutine calc_triangle_corner_at_cell_center(mesh, triangleCornerAtCellCenter)

    ! this finds the triangle corner number of a cell centre for an adjacent triangle

    type(mesh_type), intent(in) :: mesh

    integer, dimension(:,:), intent(out) :: &
         triangleCornerAtCellCenter

    integer :: &
         iCell,            & !
         iVertexOnCell,    & !
         iVertex,          & !
         iVertexDegree,    & !
         iCellVertexDegree   !

    ! loop over velocity positions
    do iCell = 1, mesh % nCells

       ! loop over triangles surrounding velocity position
       do iVertexOnCell = 1, mesh % nEdgesOnCell % array(iCell)

          ! get the vertex number associated with this cell vertex
          iVertex = mesh % verticesOnCell % array(iVertexOnCell, iCell)

          ! loop over triangle corners of this triangle
          do iVertexDegree = 1, mesh % vertexDegree

             ! get the cell number associated with this triangle corner
             iCellVertexDegree = mesh % cellsOnVertex % array(iVertexDegree, iVertex)

             ! test to see if it is the cell number of interest
             if (iCellVertexDegree == iCell) then

                triangleCornerAtCellCenter(iVertexOnCell,iCell) = iVertexDegree
                exit

             endif

          enddo ! iVertexDegree

       enddo ! iVertexOnCell

    enddo ! iCell

  end subroutine calc_triangle_corner_at_cell_center

  !-------------------------------------------------------------

  subroutine calc_edge_opposite_triangle_corner(mesh, edgeOppositeTriangleCorner, interiorVertex)

    type(mesh_type), intent(in) :: mesh
    
    integer, dimension(:,:), intent(out) :: &
         edgeOppositeTriangleCorner

    integer, dimension(:), intent(in) :: &
         interiorVertex

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell, &
         jVertexDegree, &
         jEdge, &
         jEdgeOnCell

    integer, dimension(3) :: &
         edgeMatched, &
         jEdgeVertexDegree

    ! loop over all triangles
    do iVertex = 1, mesh % nVertices

       edgeOppositeTriangleCorner(:,iVertex) = 0

       ! check are at an interior triangle
       if (interiorVertex(iVertex) == 1) then

          ! find opposite edge for each corner of triangle
          do iVertexDegree = 1, mesh % vertexDegree

             ! get cell number of this triangle corner
             iCell = mesh % cellsOnVertex % array(iVertexDegree,iVertex)

             ! initialize the matched array to none
             edgeMatched = 0

             ! check each triangle corner to see if NOT on iCell
             do jVertexDegree = 1, mesh % vertexDegree

                ! edge of this triangle corner
                jEdgeVertexDegree(jVertexDegree) = mesh % edgesOnVertex % array(jVertexDegree, iVertex)

                ! loop over edges on this cell see whether match with 
                do jEdge = 1, mesh % nEdgesOnCell % array(iCell)

                   ! edge number of this edge on the cell
                   jEdgeOnCell = mesh % edgesOnCell % array(jEdge,iCell)

                   ! see if edges are the same
                   if (jEdgeVertexDegree(jVertexDegree) == jEdgeOnCell) then

                      ! matched the vertex edge so this ISNT the right edge
                      edgeMatched(jVertexDegree) = 1

                   endif

                enddo ! iEdgeOnCell

             enddo ! jVertexDegree

             ! check have matched right number
             if (sum(edgeMatched) /= 2) then
                write(*,*) "calc_edge_opposite_triangle_corner: matched failed"
                stop
             endif

             ! find the right nonmatched edge
             do jVertexDegree = 1, mesh % vertexDegree

                if (edgeMatched(jVertexDegree) == 0) then
                   edgeOppositeTriangleCorner(iVertexDegree,iVertex) = jEdgeVertexDegree(jVertexDegree)
                endif

             enddo ! jVertexDegree

          enddo ! iVertexDegree

       endif ! interiorVertex

    enddo ! iVertex

  end subroutine calc_edge_opposite_triangle_corner

  !-------------------------------------------------------------

  subroutine basis_gradients_planar(mesh, basisGradientU, basisGradientV, interiorVertex)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         basisGradientU, &
         basisGradientV

    integer, dimension(:), intent(in) :: &
         interiorVertex

    integer :: &
         iVertex,       &
         iVertexDegree

    integer, dimension(3) :: &
         iCells

    write(*,*) "basis_gradients_planar"

    ! loop over the triangles
    do iVertex = 1, mesh % nVertices

       if (interiorVertex(iVertex) == 1) then

          ! loop over corners of triangle
          do iVertexDegree = 1, mesh % vertexDegree
             
             iCells(iVertexDegree) = mesh % cellsOnVertex % array(iVertexDegree, iVertex)
             
          enddo ! iVertexDegree
          
          call basis_gradients_planar_single(&
               basisGradientU(:,iVertex), &
               basisGradientV(:,iVertex), &
               mesh % xVertex % array(iVertex), mesh % yVertex % array(iVertex), &
               mesh % xCell % array(iCells(1)), mesh % yCell % array(iCells(1)), &
               mesh % xCell % array(iCells(2)), mesh % yCell % array(iCells(2)), &
               mesh % xCell % array(iCells(3)), mesh % yCell % array(iCells(3)))
          
       endif ! iVertex

    enddo

  end subroutine basis_gradients_planar

  !-------------------------------------------------------------

  subroutine basis_gradients_planar_single(basisGradientU, basisGradientV, &
       x0, y0, x1, y1, x2, y2, x3, y3)

    use cice_dynamics_shared, only: solve_linear_basis_system

    real(kind=RKIND), dimension(:), intent(out) :: &
         basisGradientU, &
         basisGradientV

    real(kind=RKIND), intent(in) :: &
         x0, y0, &
         x1, y1, &
         x2, y2, &
         x3, y3

    real(kind=RKIND), dimension(3,3) :: &
         leftMatrix

    real(kind=RKIND), dimension(3) :: &
         rightHandSide, &
         solutionVector

    integer :: &
         iTriangleCorner

    leftMatrix(1,1) = x1 - x0
    leftMatrix(1,2) = y1 - y0
    leftMatrix(1,3) = 1.0_RKIND

    leftMatrix(2,1) = x2 - x0
    leftMatrix(2,2) = y2 - y0
    leftMatrix(2,3) = 1.0_RKIND

    leftMatrix(3,1) = x3 - x0
    leftMatrix(3,2) = y3 - y0
    leftMatrix(3,3) = 1.0_RKIND
       
    ! loop over corners of triangle
    do iTriangleCorner = 1, 3

       rightHandSide(:)               = 0.0_RKIND
       rightHandSide(iTriangleCorner) = 1.0_RKIND

       call solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

       basisGradientU(iTriangleCorner) = solutionVector(1)
       basisGradientV(iTriangleCorner) = solutionVector(2)

    enddo ! iVertexDegree

  end subroutine basis_gradients_planar_single

  !-------------------------------------------------------------

  subroutine basis_integrals(mesh, basisIntegral, edgeOppositeTriangleCorner)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         basisIntegral

    integer, dimension(:,:), intent(in) :: &
         edgeOppositeTriangleCorner

    integer :: &
         iVertex, &
         iVertexDegree, &
         iEdge

    real(kind=RKIND) :: &
         a, b, c

    write(*,*) "basis_integrals"
    
    ! loop over triangles
    do iVertex = 1, mesh % nVertices

       do iVertexDegree = 1, mesh % vertexDegree       

          iEdge = edgeOppositeTriangleCorner(modulo(iVertexDegree-1,3)+1, iVertex)
          c = mesh % dcEdge % array(iEdge)
          
          iEdge = edgeOppositeTriangleCorner(modulo(iVertexDegree,3)+1, iVertex)
          a = mesh % dcEdge % array(iEdge)

          iEdge = edgeOppositeTriangleCorner(modulo(iVertexDegree-2,3)+1, iVertex)
          b = mesh % dcEdge % array(iEdge)          

          basisIntegral(iVertexDegree,iVertex) = &
               basis_integral(a, b, c, mesh % areaTriangle % array(iVertex))

          !write(*,*) iVertex, iVertexDegree, basisIntegral(iVertexDegree,iVertex)

       enddo ! iVertexDegree

    enddo ! iVertex

  end subroutine basis_integrals

  !-------------------------------------------------------------

  function basis_integral(a, b, c, areaTriangle) result(integral)

    real(kind=RKIND), intent(in) :: &
         a, &
         b, &
         c, &
         areaTriangle

    real(kind=RKIND) :: &
         integral

    real(kind=RKIND) :: &
         d, &
         cosGamma

    cosGamma = (a**2 + b**2 - c**2) / (2.0_RKIND * a * b)
    
    d = ((a*b) / c) * sqrt(1.0_RKIND - cosGamma**2)

    integral = (c*d) / 6.0_RKIND

  end function basis_integral

  !-------------------------------------------------------------
  ! time step
  !-------------------------------------------------------------
  
  subroutine run_dynamics_evp(mesh, icestate, dynamics, basis, forcing, boundary)

    use cice_dynamics_shared, only: ice_strength, air_stress, ocean_stress, surface_tilt, solve_velocity

    type(mesh_type), intent(in) :: mesh

    type(icestate_type), pointer :: icestate
    type(dynamics_type), pointer :: dynamics
    type(basis_type),    pointer :: basis
    type(forcing_type),  pointer :: forcing
    type(boundary_type), pointer :: boundary

    real(kind=RKIND) :: timeStep, elasticTimeStep

    call ice_strength(mesh, &     
         dynamics % icePressure % array, &
         icestate % iceAreaCell % array, &
         icestate % iceVolumeCell % array)
    
    call stress_tensor(mesh, &
         dynamics % stress1 % array,         &
         dynamics % stress2 % array,         &
         dynamics % stress12 % array,        &
         dynamics % icePressure % array,     &
         dynamics % uVelocity % array,       &
         dynamics % vVelocity % array,       &
         basis % basisGradientU % array,     &
         basis % basisGradientV % array,     &
         boundary % interiorVertex2 % array, &
         timeStep, elasticTimeStep)
    
    call divergence_stress_tensor(mesh, &
         dynamics % stressDivergenceU % array,       &
         dynamics % stressDivergenceV % array,       &
         dynamics % stress1 % array,                 &
         dynamics % stress2 % array,                 &
         dynamics % stress12 % array,                &
         basis % basisGradientU % array,             &   
         basis % basisGradientV % array,             &
         basis % basisIntegrals % array,             &
         basis % triangleCornerAtCellCenter % array, &
         boundary % boundaryCell % array,            &
         boundary % interiorVertex2 % array)

    call air_stress(&
         dynamics % airStressU % array,  &
         dynamics % airStressV % array,  &
         forcing % uAirVelocity % array, & 
         forcing % vAirVelocity % array, &
         icestate % iceAreaVertex % array)
    
    !call ocean_stress(& 
    !     dynamics % oceanStressU % array,     &
    !     dynamics % oceanStressV % array,     &
    !     dynamics % oceanStressCoeff % array, &
    !     forcing % uOceanVelocity % array,    & 
    !     forcing % vOceanVelocity % array,    &
    !     dynamics % uVelocity % array,        &   
    !     dynamics % vVelocity % array,        &
    !     icestate % iceAreaCell % array)
    
    call surface_tilt(&
         dynamics % surfaceTiltForceU % array, &
         dynamics % surfaceTiltForceV % array)
    
    call solve_velocity(mesh % nCells,         &
         boundary % interiorVertex % array,    &
         dynamics % uVelocity % array,         &   
         dynamics % vVelocity % array,         &
         iceState % totalMassCell % array,     & 
         mesh % fCell % array,                 &
         dynamics % stressDivergenceU % array, & 
         dynamics % stressDivergenceV % array, &
         dynamics % airStressU % array,        &
         dynamics % airStressV % array,        &
         dynamics % surfaceTiltForceU % array, &
         dynamics % surfaceTiltForceV % array, &
         dynamics % oceanStressU % array,      &
         dynamics % oceanStressV % array,      &
         dynamics % oceanStressCoeff % array,  &
         timeStep)

  end subroutine run_dynamics_evp

  !-------------------------------------------------------------
  
  subroutine strain_tensor(mesh, &
                           stress11,        stress12,      &
                           stress12,                       &
                           icePressure,                    &
                           uVelocity,      vVelocity,      &
                           basisGradientU, basisGradientV, &
                           interiorVertex,                 &
                           timeStep,       elasticTimeStep)

    use cice_dynamics_shared, only: evp_constitutive_relation

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         stress1, &
         stress2, &
         stress12

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         basisGradientU, &
         basisGradientV

    real(kind=RKIND), dimension(:), intent(in) :: &
         icePressure

    integer, dimension(:), intent(in) :: &
         interiorVertex

    real(kind=RKIND), intent(in) :: &
         timeStep,  &
         elasticTimeStep

    real(kind=RKIND) :: &
         strainDivergence, & ! divergence of the strain rate tensor
         strainTension,    & ! horizontal tension of the strain rate tensor
         strainShearing      ! horizontal shearing of the strain rate tensor

    integer :: &
         iVertex,       &
         iVertexDegree, &
         iCell

    ! loop over the triangles
    do iVertex = 1, mesh % nVertices

       if (interiorVertex(iVertex) == 1) then
          
          strainDivergence = 0.0_RKIND
          strainTension    = 0.0_RKIND
          strainShearing   = 0.0_RKIND
          
          ! loop over corners of triangle
          do iVertexDegree = 1, mesh % vertexDegree
             
             iCell = mesh % cellsOnVertex % array(iVertexDegree, iVertex)
             
             ! divergence of the strain rate tensor
             strainDivergence = strainDivergence + &
                  uVelocity(iCell) * basisGradientU(iVertexDegree,iVertex) + &
                  vVelocity(iCell) * basisGradientV(iVertexDegree,iVertex)
             
             ! horizontal tension of the strain rate tensor
             strainTension = strainTension + &
                  uVelocity(iCell) * basisGradientU(iVertexDegree,iVertex) - &
                  vVelocity(iCell) * basisGradientV(iVertexDegree,iVertex)
             
             ! horizontal shearing of the strain rate tensor
             strainShearing = strainShearing + &
                  uVelocity(iCell) * basisGradientV(iVertexDegree,iVertex) + &
                  vVelocity(iCell) * basisGradientU(iVertexDegree,iVertex)
             
          enddo ! iVertexDegree
          
          ! loop over corners of triangle
          do iVertexDegree = 1, mesh % vertexDegree
             
             iCell = mesh % cellsOnVertex % array(iVertexDegree, iVertex)
             
             !call evp_constitutive_relation(stress1(iVertexDegree,iVertex),                     &
             !                               stress2(iVertexDegree,iVertex),                     &
             !                               stress12(iVertexDegree,iVertex),                    &
             !                               strainDivergence, strainTension, strainShearing,    &
             !                               icePressure(iCell), mesh % areaCell % array(iCell), &
             !                               timeStep, elasticTimeStep)
  
          enddo ! iVertexDegree

       endif ! interiorVertex

    enddo ! iVertex

  end subroutine stress_tensor

  !-------------------------------------------------------------

  subroutine divergence_stress_tensor(mesh, &
                                      divergenceStressU, divergenceStressV, &
                                      stress1,            stress2,          &
                                      stress12,                             &
                                      basisGradientU,     basisGradientV,   &
                                      basisIntegrals,                       &
                                      triangleCornerAtCellCenter,           &
                                      boundaryCell,       interiorVertex2)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         divergenceStressU, &
         divergenceStressV

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         stress1,        &
         stress2,        &
         stress12,       &
         basisGradientU, &
         basisGradientV, &
         basisIntegrals

    integer, dimension(:,:), intent(in) :: &
         triangleCornerAtCellCenter

    integer, dimension(:), intent(in) :: &
         boundaryCell, &
         interiorVertex2

    integer :: &
         iCell,         &
         iVertexOnCell, &
         iVertex,       &
         iVertexDegree, &
         jVertexDegree


    real(kind=RKIND) :: &
         stress1contribU, &
         stress2contribU, &
         stress12contribU, &
         stress1contribV, &
         stress2contribV, &
         stress12contribV

    write(*,*) "divergence start", mesh % nCells

    ! loop over velocity positions
    do iCell = 1, mesh % nCells

       ! only cells interior to our smaller domain
       if (boundaryCell(iCell) == 0) then

          divergenceStressU(iCell) = 0.0_RKIND
          divergenceStressV(iCell) = 0.0_RKIND   
          
          ! loop over triangles surrounding velocity position
          do iVertexOnCell = 1, mesh % nEdgesOnCell % array(iCell)
             
             ! vertex index of the current triangle
             iVertex = mesh % verticesOnCell % array(iVertexOnCell, iCell)

             stress1contribU  = 0.0_RKIND
             stress2contribU  = 0.0_RKIND
             stress12contribU = 0.0_RKIND
             
             stress1contribV  = 0.0_RKIND
             stress2contribV  = 0.0_RKIND
             stress12contribV = 0.0_RKIND
             
             ! only internal triangles at the mo
             !if (interiorVertex2(iVertex) == 1) then
             
                ! vertex number of vertex triangle from this cell
                jVertexDegree = triangleCornerAtCellCenter(iVertexOnCell, iCell)
                
                do iVertexDegree = 1, mesh % vertexDegree

                   stress1contribU  = stress1contribU  + stress1(iVertexDegree,iVertex)  * basisGradientU(iVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex)
                   stress2contribU  = 0.0_RKIND
                   stress12contribU = stress12contribU + stress12(iVertexDegree,iVertex) * basisGradientV(iVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex)
                   
                   stress1contribV  = 0.0_RKIND
                   stress2contribV  = stress2contribV  + stress2(iVertexDegree,iVertex)  * basisGradientV(jVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex)
                   stress12contribV = stress12contribV + stress12(iVertexDegree,iVertex) * basisGradientU(jVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex)
                
                   ! divergence of the stress tensor in u direction
                   !divergenceStressU(iCell) = divergenceStressU(iCell) + &
                   !     stress1(iVertexDegree,iVertex)  * basisGradientU(iVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex) * 0.5_RKIND + &
                   !     stress2(iVertexDegree,iVertex)  * basisGradientU(iVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex) * 0.5_RKIND + &
                   !     stress12(iVertexDegree,iVertex) * basisGradientV(iVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex)
                   
                   ! divergence of the stress tensor in v direction
                   !divergenceStressV(iCell) = divergenceStressV(iCell) + &
                   !     stress1(iVertexDegree,iVertex)  * basisGradientV(jVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex) * 0.5_RKIND + &
                   !     stress2(iVertexDegree,iVertex)  * basisGradientV(jVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex) * 0.5_RKIND + &
                   !     stress12(iVertexDegree,iVertex) * basisGradientU(jVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex)

                   write(*,*) iCell, iVertexOnCell, iVertexDegree, stress1(iVertexDegree,iVertex), basisGradientU(iVertexDegree,iVertex), basisIntegrals(iVertexDegree,iVertex)
                   
                enddo ! iVertexDegree
                
             !endif ! interiorVertex2

             divergenceStressU(iCell) = divergenceStressU(iCell) + &
                  (stress1contribU + stress12contribU) / mesh % areaTriangle % array(iVertex)

             divergenceStressV(iCell) = divergenceStressV(iCell) + &
                  (stress2contribU + stress12contribU) / mesh % areaTriangle % array(iVertex)

          enddo ! iVertexOnCell

          write(*,*) iCell, divergenceStressU(iCell), divergenceStressV(iCell)

       endif ! boundaryCell

    enddo ! iCell

    write(*,*) "divergence"

  end subroutine divergence_stress_tensor

  !-------------------------------------------------------------

  subroutine test_tri_divergence_stress_tensor(mesh, dynamics, basis, boundary)

    use cice_testing, only: divergence_stress_test_stress_set_tri

    type(mesh_type), intent(in) :: mesh
    type(dynamics_type), pointer :: dynamics
    type(basis_type),    pointer :: basis
    type(boundary_type), pointer :: boundary

    write(*,*) "divergence test"

    call divergence_stress_test_stress_set_tri(mesh, &
                                               dynamics % stress1 % array, &
                                               dynamics % stress2 % array, &
                                               dynamics % stress12 % array)

    call divergence_stress_tensor(mesh, &
         dynamics % stressDivergenceU % array,       &
         dynamics % stressDivergenceV % array,       &
         dynamics % stress1 % array,                 &
         dynamics % stress2 % array,                 &
         dynamics % stress12 % array,                &
         basis % basisGradientU % array,             &
         basis % basisGradientV % array,             &
         basis % basisIntegrals % array,             &
         basis % triangleCornerAtCellCenter % array, &
         boundary % boundaryCell % array,            &
         boundary % interiorVertex % array)

    stop

  end subroutine test_tri_divergence_stress_tensor

  !-------------------------------------------------------------
  ! boundary
  !-------------------------------------------------------------

  subroutine divergence_stress_tensor_boundary(mesh, &
                                               divergenceStressU, divergenceStressV, &
                                               uVelocity,         vVelocity,         &
                                               icePressure,       basisIntegrals,    &
                                               basisGradientU,    basisGradientV,    &
                                               boundaryCell,      boundaryCell2,     &
                                               boundaryVertex2,                      &
                                               timeStep,          elasticTimeStep)

    use cice_dynamics_shared, only: evp_constitutive_relation

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(inout) :: &
         divergenceStressU, &
         divergenceStressV

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, &
         vVelocity, &
         icePressure

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         basisIntegrals, &
         basisGradientU, &
         basisGradientV

    integer, dimension(:), intent(in) :: &
          boundaryCell,    &
          boundaryCell2,   &
          boundaryVertex2

    real(kind=RKIND), intent(in) :: &
         timeStep,  &
         elasticTimeStep

    real(kind=RKIND) :: &
         strainDivergence, &
         strainTension, &
         strainShearing, &
         stress1, &
         stress2, &
         stress12, &
         divergenceStressUBoundary, &
         divergenceStressVBoundary, &
         areaCorrection

    integer :: &
         iCell, &
         iVertexOnCell, &
         iVertex, &
         iVertexDegree, &
         iCellTriangleCorner, &
         jVertexDegree

    ! loop over all velocity points
    do iCell = 1, mesh % nCells

       ! interior boundary velocity points
       if (boundaryCell2(iCell) == 1) then

          ! find the triangles we still need to sum
          do iVertexOnCell = 1, mesh % nEdgesOnCell % array(iCell) 

             iVertex = mesh % verticesOnCell % array(iVertexOnCell,iCell)

             if (boundaryVertex2(iVertex) == 1) then

                strainDivergence = 0.0_RKIND
                strainTension    = 0.0_RKIND
                strainShearing   = 0.0_RKIND

                areaCorrection = 0.0_RKIND

                ! loop over triangle vertices
                do iVertexDegree = 1, mesh % vertexDegree

                   iCellTriangleCorner = mesh % cellsOnVertex % array(iVertexDegree,iVertex)

                   if (boundaryCell(iCellTriangleCorner) == 0) then
                      ! interior triangle corner
                      areaCorrection = areaCorrection + mesh % kiteAreasOnVertex % array(iVertexDegree,iVertex)

                      ! divergence of the strain rate tensor
                      strainDivergence = strainDivergence + &
                           uVelocity(iCellTriangleCorner) * basisGradientU(iVertexDegree,iVertex) + &
                           vVelocity(iCellTriangleCorner) * basisGradientV(iVertexDegree,iVertex)
                      
                      ! horizontal tension of the strain rate tensor
                      strainTension = strainTension + &
                           uVelocity(iCellTriangleCorner) * basisGradientU(iVertexDegree,iVertex) - &
                           vVelocity(iCellTriangleCorner) * basisGradientV(iVertexDegree,iVertex)
                      
                      ! horizontal shearing of the strain rate tensor
                      strainShearing = strainShearing + &
                           uVelocity(iCellTriangleCorner) * basisGradientV(iVertexDegree,iVertex) + &
                           vVelocity(iCellTriangleCorner) * basisGradientU(iVertexDegree,iVertex)

                   else
                      ! exterior triangle corner

                      ! divergence of the strain rate tensor
                      strainDivergence = strainDivergence - &
                           uVelocity(iCell) * basisGradientU(iVertexDegree,iVertex) - &
                           vVelocity(iCell) * basisGradientV(iVertexDegree,iVertex)
                      
                      ! horizontal tension of the strain rate tensor
                      strainTension = strainTension - &
                           uVelocity(iCell) * basisGradientU(iVertexDegree,iVertex) + &
                           vVelocity(iCell) * basisGradientV(iVertexDegree,iVertex)
                      
                      ! horizontal shearing of the strain rate tensor
                      strainShearing = strainShearing - &
                           uVelocity(iCell) * basisGradientV(iVertexDegree,iVertex) - &
                           vVelocity(iCell) * basisGradientU(iVertexDegree,iVertex)

                   endif ! boundaryCell
                   
                enddo ! iVertexDegree

                areaCorrection = areaCorrection / mesh % areaCell % array(iCell)
                
                ! calculate stress assuming all ice strength equal that at iCell
                !call evp_constitutive_relation(stress1,          stress2,       stress12,          &
                !                               strainDivergence, strainTension, strainShearing,    &
                !                               icePressure(iCell), mesh % areaCell % array(iCell), &
                !                               timeStep, elasticTimeStep)

                ! vertex number of vertex triangle from this cell
                !jVertexDegree = triangleCornerAtCellCenter(iVertexOnCell, iCell) ????
                
                divergenceStressUBoundary = 0.0_RKIND
                divergenceStressVBoundary = 0.0_RKIND

                do iVertexDegree = 1, mesh % vertexDegree
                   
                   ! divergence of the stress tensor in u direction
                   divergenceStressUBoundary = divergenceStressUBoundary + &
                        stress1  * basisGradientU(iVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex) * 0.5_RKIND + &
                        stress2  * basisGradientU(iVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex) * 0.5_RKIND + &
                        stress12 * basisGradientV(iVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex)
                   
                   ! divergence of the stress tensor in v direction
                   divergenceStressVBoundary = divergenceStressVBoundary + &
                        stress1  * basisGradientV(jVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex) * 0.5_RKIND + &
                        stress2  * basisGradientV(jVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex) * 0.5_RKIND + &
                        stress12 * basisGradientU(jVertexDegree,iVertex) * basisIntegrals(iVertexDegree,iVertex)
                   
                enddo ! iVertexDegree

                divergenceStressU(iCell) = divergenceStressU(iCell) + divergenceStressUBoundary * areaCorrection
                divergenceStressV(iCell) = divergenceStressV(iCell) + divergenceStressVBoundary * areaCorrection

             endif ! boundaryVertex2

          enddo ! iVertexOnCell

       endif ! boundaryCell2
       
    enddo ! iCell

  end subroutine divergence_stress_tensor_boundary

  !------------------------------------------------------------- 

  subroutine set_boundary_velocity(mesh, &
                                   uVelocityBoundary, vVelocityBoundary, &
                                   uVelocity,         vVelocity,         &
                                   boundaryCell)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         uVelocityBoundary, &
         vVelocityBoundary

    real(kind=RKIND), dimension(:), intent(in) :: &    
         uVelocity, &
         vVelocity

    integer, dimension(:), intent(in) :: &
         boundaryCell

    integer :: &
         iCell, &
         iEdgeOnCell, &
         iCellNeighbour

    ! loop over all possible velocity points
    do iCell = 1, mesh % nCells

       ! only boundary velocity points
       if (boundaryCell(iCell) == 1) then

          ! loop over surrounding cells
          do iEdgeOnCell = 1, mesh % nEdgesOnCell % array(iCell)

             ! cell number of neighbouring cell
             iCellNeighbour = mesh % cellsOnCell % array(iEdgeOnCell,iCell)

             ! check if this is interior cell
             if (boundaryCell(iCellNeighbour) == 0) then

                uVelocityBoundary(iEdgeOnCell,iCell) = -uVelocity(iCellNeighbour)
                vVelocityBoundary(iEdgeOnCell,iCell) = -vVelocity(iCellNeighbour)

             endif ! boundaryCell

          enddo ! iEdgeOnCell

       endif ! boundaryCell

    enddo ! iCell

  end subroutine set_boundary_velocity

  !-------------------------------------------------------------

end module cice_dynamics_evp_tri

