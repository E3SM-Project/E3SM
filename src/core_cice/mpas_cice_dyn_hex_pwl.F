module cice_dyn_hex_pwl

  use mpas_grid_types

  implicit none

  private
  public :: &
       init_dynamics_hex_pwl, &
       internal_stress_pwl, &
       strain_tensor_hex_pwl, &
       stress_divergence_hex_pwl

contains

  !-------------------------------------------------------------
  ! Initialization
  !-------------------------------------------------------------

  subroutine init_dynamics_hex_pwl(mesh, hexvar, hexpwl, boundary)

    use cice_mesh, only: calc_cell_vertices_at_vertex

    type(MPAS_pool_type), intent(in) :: mesh

    type(MPAS_pool_type), pointer :: hexvar
    type(MPAS_pool_type), pointer :: boundary
    type(MPAS_pool_type), pointer :: hexpwl

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex

    real(kind=RKIND), dimension(:,:), pointer :: &
         basisSubArea

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisGradientU, &
         basisGradientV

    call MPAS_pool_get_array(hexvar, "cellVerticesAtVertex", cellVerticesAtVertex)
    call MPAS_pool_get_array(hexpwl, "basisSubArea", basisSubArea)
    call MPAS_pool_get_array(hexpwl, "basisGradientU", basisGradientU)
    call MPAS_pool_get_array(hexpwl, "basisGradientV", basisGradientV)

    call calc_cell_vertices_at_vertex(mesh, &
         cellVerticesAtVertex)

    call init_dynamics_hex_pwl_basis(mesh, &
         basisSubArea, &
         basisGradientU, &
         basisGradientV)

  end subroutine init_dynamics_hex_pwl

  !-------------------------------------------------------------

  subroutine init_dynamics_hex_pwl_basis(mesh,           basisSubArea, &
                                         basisGradientU, basisGradientV)

    use cice_numerics, only: solve_linear_basis_system

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         basisSubArea
         
    real(kind=RKIND), dimension(:,:,:), intent(out) :: &    
         basisGradientU, &
         basisGradientV

    real(kind=RKIND) :: &
         xPWLCentre, &
         yPWLCentre, &
         alphaPWL, &
         a, b, c, s

    integer :: &
         iCell, &
         iVertexOnCell, &
         iVertex, &
         iEdgeOnCell, &
         iEdge, &
         iVertexOnCell1, &
         iVertexOnCell2, &
         iVertex1, &
         iVertex2

    real(kind=RKIND), dimension(3,3) :: &
         leftMatrix

    real(kind=RKIND), dimension(3) :: &
         rightHandSide, &
         solutionVector

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell, &
         edgesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         dvEdge

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "dvEdge", dvEdge)

    ! loop over cells
    do iCell = 1, nCells
       
       alphaPWL = 1.0_RKIND / real(nEdgesOnCell(iCell),RKIND)

       ! determine cell centre for piecewise linear basis
       xPWLCentre = 0.0_RKIND
       yPWLCentre = 0.0_RKIND

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          iVertex = verticesOnCell(iVertexOnCell,iCell)

          xPWLCentre = xPWLCentre + alphaPWL * xVertex(iVertex)
          yPWLCentre = yPWLCentre + alphaPWL * yVertex(iVertex)

       enddo ! iVertexOnCell

       ! calculate the area of the subcells
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          iEdge = edgesOnCell(iEdgeOnCell,iCell)
          iVertexOnCell1 = iEdgeOnCell
          iVertexOnCell2 = wrapped_index(iEdgeOnCell + 1, nEdgesOnCell(iCell))

          iVertex1 = verticesOnCell(iVertexOnCell1,iCell)
          iVertex2 = verticesOnCell(iVertexOnCell2,iCell)

          c = dvEdge(iEdge)
          a = sqrt((xVertex(iVertex1) - xPWLCentre)**2 + &
                   (yVertex(iVertex1) - yPWLCentre)**2)
          b = sqrt((xVertex(iVertex2) - xPWLCentre)**2 + &
                   (yVertex(iVertex2) - yPWLCentre)**2)

          s = (a + b + c) * 0.5_RKIND

          ! Heron's formula 
          basisSubArea(iEdgeOnCell,iCell) = &
               sqrt(s * (s-a) * (s-b) * (s-c))

          ! set up left hand matrix
          leftMatrix(1,1) = xVertex(iVertex1) - xPWLCentre
          leftMatrix(1,2) = yVertex(iVertex1) - yPWLCentre
          leftMatrix(1,3) = 1.0_RKIND
          
          leftMatrix(2,1) = xVertex(iVertex2) - xPWLCentre
          leftMatrix(2,2) = yVertex(iVertex2) - yPWLCentre
          leftMatrix(2,3) = 1.0_RKIND
          
          leftMatrix(3,1) = 0.0_RKIND
          leftMatrix(3,2) = 0.0_RKIND
          leftMatrix(3,3) = 1.0_RKIND

          ! first basis
          rightHandSide(1) = 1.0_RKIND
          rightHandSide(2) = 0.0_RKIND
          rightHandSide(3) = 0.0_RKIND

          call solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

          basisGradientU(1,iEdgeOnCell,iCell) = solutionVector(1)
          basisGradientV(1,iEdgeOnCell,iCell) = solutionVector(2)

          ! second basis
          rightHandSide(1) = 0.0_RKIND
          rightHandSide(2) = 1.0_RKIND
          rightHandSide(3) = 0.0_RKIND

          call solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

          basisGradientU(2,iEdgeOnCell,iCell) = solutionVector(1)
          basisGradientV(2,iEdgeOnCell,iCell) = solutionVector(2)

       enddo ! iEdgeOnCell

    enddo ! iCell

  end subroutine init_dynamics_hex_pwl_basis

  !-------------------------------------------------------------
  ! time step
  !-------------------------------------------------------------

  subroutine internal_stress_pwl(&
       mesh, &
       hexvar,              hexpwl, &
       uVelocity,           vVelocity, &
       icePressure,         replacementPressure, &
       dtElastic,           evpDamping, &
       stressDivergenceU,   stressDivergenceV, &
       solveStress,         solveVelocity)

    use cice_dyn_hex_variational, only: &
         stress_tensor_hex_variational

    type (MPAS_pool_type), pointer :: &
         mesh, &
         hexvar, &
         hexpwl

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, &
         vVelocity       

    real(kind=RKIND), dimension(:), intent(inout) :: &
         icePressure, &
         replacementPressure

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

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11, &
         strain22, &
         strain12, &
         stress11, &
         stress22, &
         stress12, &
         basisSubArea

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisGradientUPWL, &
         basisGradientVPWL

    call MPAS_pool_get_array(hexvar, "strain11", strain11)
    call MPAS_pool_get_array(hexvar, "strain22", strain22)
    call MPAS_pool_get_array(hexvar, "strain12", strain12)
    call MPAS_pool_get_array(hexvar, "stress11", stress11)
    call MPAS_pool_get_array(hexvar, "stress22", stress22)
    call MPAS_pool_get_array(hexvar, "stress12", stress12)
    call MPAS_pool_get_array(hexvar, "cellVerticesAtVertex", cellVerticesAtVertex)

    call MPAS_pool_get_array(hexpwl, "basisGradientU", basisGradientUPWL)
    call MPAS_pool_get_array(hexpwl, "basisGradientV", basisGradientVPWL)
    call MPAS_pool_get_array(hexpwl, "basisSubArea", basisSubArea)

    call strain_tensor_hex_pwl(mesh, &
         strain11,          &
         strain22,          &
         strain12,          &
         uVelocity,         &
         vVelocity,         &
         basisGradientUPWL, &
         basisGradientVPWL, &
         solveStress)
    
    call stress_tensor_hex_variational(mesh, &
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
    
    call stress_divergence_hex_pwl(mesh, &
         stressDivergenceU,    &
         stressDivergenceV,    &
         stress11,             &
         stress22,             &
         stress12,             &
         basisGradientUPWL,    &
         basisGradientVPWL,    &
         basisSubArea,         &
         cellVerticesAtVertex, &
         solveVelocity)

  end subroutine internal_stress_pwl

  !-------------------------------------------------------------
  
  subroutine strain_tensor_hex_pwl(mesh, &
                                   strain11,       strain22,       &
                                   strain12,                       &
                                   uVelocity,      vVelocity,      &
                                   basisGradientU, basisGradientV, &
                                   solveStress)

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         strain11, &
         strain22, &
         strain12
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         basisGradientU, &
         basisGradientV
    
    integer, dimension(:), intent(in) :: &
         solveStress

    real(kind=RKIND) :: &
         alphaPWL

    integer :: &
         iCell, &
         jVertexOnCell, &
         iVertexOnCell, &
         iVertex, &
         jVertex, &
         iSubCell1, &
         iSubCell2

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)

    ! loop over cells
    do iCell = 1, nCells

       strain11(:,iCell) = 0.0_RKIND
       strain22(:,iCell) = 0.0_RKIND
       strain12(:,iCell) = 0.0_RKIND

       if (solveStress(iCell) == 1) then

          alphaPWL = 1.0_RKIND / real(nEdgesOnCell(iCell), RKIND)

          ! loop over velocity points surrounding cell - location of stress and derivative
          do jVertexOnCell = 1, nEdgesOnCell(iCell)

             jVertex = verticesOnCell(jVertexOnCell,iCell)

             iSubCell1 = jVertexOnCell
             iSubCell2 = wrapped_index(jVertexOnCell + 1, nEdgesOnCell(iCell))

             ! loop over basis functions
             do iVertexOnCell = 1, nEdgesOnCell(iCell)

                iVertex = verticesOnCell(iVertexOnCell,iCell)

                strain11(jVertexOnCell,iCell) = strain11(jVertexOnCell,iCell) + &
                     uVelocity(iVertex) * (pwl_basis_gradient(nEdgesOnCell(iCell), basisGradientU(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell1)  + &
                                           pwl_basis_gradient(nEdgesOnCell(iCell), basisGradientU(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell2))

                strain22(jVertexOnCell,iCell) = strain22(jVertexOnCell,iCell) + &
                     vVelocity(iVertex) * (pwl_basis_gradient(nEdgesOnCell(iCell), basisGradientV(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell1)  + &
                                           pwl_basis_gradient(nEdgesOnCell(iCell), basisGradientV(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell2))

                strain12(jVertexOnCell,iCell) = strain12(jVertexOnCell,iCell) + 0.5_RKIND * ( &
                     uVelocity(iVertex) * (pwl_basis_gradient(nEdgesOnCell(iCell), basisGradientV(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell1)  + &
                                           pwl_basis_gradient(nEdgesOnCell(iCell), basisGradientV(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell2)) + &
                     vVelocity(iVertex) * (pwl_basis_gradient(nEdgesOnCell(iCell), basisGradientU(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell1)  + &
                                           pwl_basis_gradient(nEdgesOnCell(iCell), basisGradientU(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell2)))

             enddo ! iVertexOnCell

          enddo ! jVertexOnCell

       endif ! solveStress

    enddo ! iCell

  end subroutine strain_tensor_hex_pwl

  !-------------------------------------------------------------

  subroutine stress_divergence_hex_pwl(mesh, &
                                       stressDivergenceU,    stressDivergenceV, &
                                       stress11,             stress22,          &
                                       stress12,                                &
                                       basisGradientU,       basisGradientV,    &
                                       basisSubArea, &
                                       cellVerticesAtVertex, solveVelocity)

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, &
         stressDivergenceV

    real(kind=RKIND), dimension(:,:), intent(in) :: &    
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &   
         basisGradientU, &
         basisGradientV

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         basisSubArea

    integer, dimension(:,:), intent(in) :: &
         cellVerticesAtVertex

    integer, dimension(:), intent(in) :: &
         solveVelocity

    real(kind=RKIND) :: &
         stressDivergenceUCell, &
         stressDivergenceVCell, &
         alphaPWL

    integer :: &
         iVertex, &
         iCellOnVertex, &
         iCell, &
         iVertexOnCell, &
         jVertexOnCell, &
         iSubCell

    integer, pointer :: &
         nVertices, &
         vertexDegree

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    real(kind=RKIND), dimension(:), pointer :: &
         areaTriangle

    ! init variables
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "areaTriangle", areaTriangle)

    ! loop over velocity positions
    do iVertex = 1, nVertices

       stressDivergenceU(iVertex) = 0.0_RKIND
       stressDivergenceV(iVertex) = 0.0_RKIND

       if (solveVelocity(iVertex) == 1) then

          ! loop over surrounding cells
          do iCellOnVertex = 1, vertexDegree
             
             ! get the cell number of this cell
             iCell = cellsOnVertex(iCellOnVertex, iVertex)

             alphaPWL = 1.0_RKIND / real(nEdgesOnCell(iCell), RKIND)
             
             ! get the vertexOnCell number of the iVertex velocity point from cell iCell
             jVertexOnCell = cellVerticesAtVertex(iCellOnVertex,iVertex)

             stressDivergenceUCell = 0.0_RKIND
             stressDivergenceVCell = 0.0_RKIND

             ! loop over the vertices of the surrounding cell
             do iVertexOnCell = 1, nEdgesOnCell(iCell)

                ! loop over subcells
                do iSubCell = 1, nEdgesOnCell(iCell)

                   stressDivergenceUCell = stressDivergenceUCell + &
                        stress11(iVertexOnCell,iCell) * basisSubArea(iSubCell,iCell) * &
                                                        pwl_basis_gradient(nEdgesOnCell(iCell), basisGradientU(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell) + &
                        stress12(iVertexOnCell,iCell) * basisSubArea(iSubCell,iCell) * &
                                                        pwl_basis_gradient(nEdgesOnCell(iCell), basisGradientV(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell)
                   
                   stressDivergenceVCell = stressDivergenceVCell + &
                        stress22(iVertexOnCell,iCell) * basisSubArea(iSubCell,iCell) * &
                                                        pwl_basis_gradient(nEdgesOnCell(iCell), basisGradientV(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell) + &
                        stress12(iVertexOnCell,iCell) * basisSubArea(iSubCell,iCell) * &
                                                        pwl_basis_gradient(nEdgesOnCell(iCell), basisGradientU(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell)

                enddo ! iSubCell
                
             enddo ! jVertex

             stressDivergenceU(iVertex) = stressDivergenceU(iVertex) - stressDivergenceUCell
             stressDivergenceV(iVertex) = stressDivergenceV(iVertex) - stressDivergenceVCell
             
          enddo ! iCellOnVertex

          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) / areaTriangle(iVertex)
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) / areaTriangle(iVertex)

       endif ! solveVelocity

    enddo ! iVertex
 
  end subroutine stress_divergence_hex_pwl

  !-------------------------------------------------------------

  function pwl_basis_gradient(nEdgesOnCell, basisGradient, alphaPWL, basisVertexOnCell, iSubCell) result(grad)
    
    integer, intent(in) :: &
         nEdgesOnCell, &
         basisVertexOnCell, & ! basis function vertex
         iSubCell             ! subcell to calculate

    real(kind=RKIND), intent(in) :: &
         alphaPWL

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         basisGradient

    real(kind=RKIND) :: grad
    
    grad = basisGradient(1,iSubCell) * &
           (alphaPWL + merge(1.0_RKIND, 0.0_RKIND, iSubCell == basisVertexOnCell)) + &
           basisGradient(2,iSubCell) * &
           (alphaPWL + merge(1.0_RKIND, 0.0_RKIND, iSubCell == wrapped_index(basisVertexOnCell - 1, nEdgesOnCell)))

  end function pwl_basis_gradient

  !-------------------------------------------------------------

  function wrapped_index(input,nelements) result(output)

    integer, intent(in) :: &
         input, &
         nelements

    integer :: output

    output = mod(input - 1, nelements) + 1

  end function wrapped_index

  !-------------------------------------------------------------

end module cice_dyn_hex_pwl
