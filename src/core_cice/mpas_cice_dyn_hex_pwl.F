module cice_dyn_hex_pwl

  use mpas_grid_types

  implicit none

  private
  public :: &
       init_dynamics_hex_pwl, &
       strain_tensor_hex_pwl, &
       stress_divergence_hex_pwl

contains

  !-------------------------------------------------------------
  ! Initialization
  !-------------------------------------------------------------

  subroutine init_dynamics_hex_pwl(mesh, hexvar, hexpwl, boundary)

    use cice_mesh, only: calc_cell_vertices_at_vertex

    type(mesh_type), intent(in) :: mesh

    type(hexvar_type),   pointer :: hexvar
    type(boundary_type), pointer :: boundary
    type(hexpwl_type),   pointer :: hexpwl

    call calc_cell_vertices_at_vertex(mesh, &
         hexvar % cellVerticesAtVertex % array)

    call init_dynamics_hex_pwl_basis(mesh, &
         hexpwl % basisSubArea % array, &
         hexpwl % basisGradientU % array, &
         hexpwl % basisGradientV % array)

  end subroutine init_dynamics_hex_pwl

  !-------------------------------------------------------------

  subroutine init_dynamics_hex_pwl_basis(mesh,           basisSubArea, &
                                         basisGradientU, basisGradientV)

    use cice_numerics, only: solve_linear_basis_system

    type(mesh_type), intent(in) :: mesh

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
         nEdgesOnCell, &
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

    ! loop over cells
    do iCell = 1, mesh % nCells

       nEdgesOnCell = mesh % nEdgesOnCell % array(iCell)
       
       alphaPWL = 1.0_RKIND / real(nEdgesOnCell,RKIND)

       ! determine cell centre for piecewise linear basis
       xPWLCentre = 0.0_RKIND
       yPWLCentre = 0.0_RKIND

       do iVertexOnCell = 1, nEdgesOnCell

          iVertex = mesh % verticesOnCell % array(iVertexOnCell,iCell)

          xPWLCentre = xPWLCentre + alphaPWL * mesh % xVertex % array(iVertex)
          yPWLCentre = yPWLCentre + alphaPWL * mesh % yVertex % array(iVertex)

       enddo ! iVertexOnCell

       ! calculate the area of the subcells
       do iEdgeOnCell = 1, nEdgesOnCell

          iEdge = mesh % edgesOnCell % array(iEdgeOnCell,iCell)
          iVertexOnCell1 = iEdgeOnCell
          iVertexOnCell2 = wrapped_index(iEdgeOnCell + 1, nEdgesOnCell)

          iVertex1 = mesh % verticesOnCell % array(iVertexOnCell1,iCell)
          iVertex2 = mesh % verticesOnCell % array(iVertexOnCell2,iCell)

          c = mesh % dvEdge % array(iEdge)
          a = sqrt((mesh % xVertex % array(iVertex1) - xPWLCentre)**2 + &
                   (mesh % yVertex % array(iVertex1) - yPWLCentre)**2)
          b = sqrt((mesh % xVertex % array(iVertex2) - xPWLCentre)**2 + &
                   (mesh % yVertex % array(iVertex2) - yPWLCentre)**2)

          s = (a + b + c) * 0.5_RKIND

          ! Heron's formula 
          basisSubArea(iEdgeOnCell,iCell) = &
               sqrt(s * (s-a) * (s-b) * (s-c))

          ! set up left hand matrix
          leftMatrix(1,1) = mesh % xVertex % array(iVertex1) - xPWLCentre
          leftMatrix(1,2) = mesh % yVertex % array(iVertex1) - yPWLCentre
          leftMatrix(1,3) = 1.0_RKIND
          
          leftMatrix(2,1) = mesh % xVertex % array(iVertex2) - xPWLCentre
          leftMatrix(2,2) = mesh % yVertex % array(iVertex2) - yPWLCentre
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
  
  subroutine strain_tensor_hex_pwl(mesh, &
                                   strain11,       strain22,       &
                                   strain12,                       &
                                   uVelocity,      vVelocity,      &
                                   basisGradientU, basisGradientV, &
                                   solveStress)

    type(mesh_type), intent(in) :: mesh

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
         nEdgesOnCell, &
         jVertexOnCell, &
         iVertexOnCell, &
         iVertex, &
         jVertex, &
         iSubCell1, &
         iSubCell2

    ! loop over cells
    do iCell = 1, mesh % nCells

       strain11(:,iCell) = 0.0_RKIND
       strain22(:,iCell) = 0.0_RKIND
       strain12(:,iCell) = 0.0_RKIND

       if (solveStress(iCell) == 1) then

          nEdgesOnCell = mesh % nEdgesOnCell % array(iCell)

          alphaPWL = 1.0_RKIND / real(nEdgesOnCell, RKIND)

          ! loop over velocity points surrounding cell - location of stress and derivative
          do jVertexOnCell = 1, nEdgesOnCell

             jVertex = mesh % verticesOnCell % array(jVertexOnCell,iCell)

             iSubCell1 = jVertexOnCell
             iSubCell2 = wrapped_index(jVertexOnCell + 1, nEdgesOnCell)

             ! loop over basis functions
             do iVertexOnCell = 1, nEdgesOnCell

                iVertex = mesh % verticesOnCell % array(iVertexOnCell,iCell)

                strain11(jVertexOnCell,iCell) = strain11(jVertexOnCell,iCell) + &
                     uVelocity(iVertex) * (pwl_basis_gradient(nEdgesOnCell, basisGradientU(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell1)  + &
                                           pwl_basis_gradient(nEdgesOnCell, basisGradientU(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell2))

                strain22(jVertexOnCell,iCell) = strain22(jVertexOnCell,iCell) + &
                     vVelocity(iVertex) * (pwl_basis_gradient(nEdgesOnCell, basisGradientV(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell1)  + &
                                           pwl_basis_gradient(nEdgesOnCell, basisGradientV(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell2))

                strain12(jVertexOnCell,iCell) = strain12(jVertexOnCell,iCell) + 0.5_RKIND * ( &
                     uVelocity(iVertex) * (pwl_basis_gradient(nEdgesOnCell, basisGradientV(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell1)  + &
                                           pwl_basis_gradient(nEdgesOnCell, basisGradientV(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell2)) + &
                     vVelocity(iVertex) * (pwl_basis_gradient(nEdgesOnCell, basisGradientU(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell1)  + &
                                           pwl_basis_gradient(nEdgesOnCell, basisGradientU(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell2)))

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

    type(mesh_type), intent(in) :: mesh

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
         nEdgesOnCell, &
         iVertexOnCell, &
         jVertexOnCell, &
         iSubCell

    ! loop over velocity positions
    do iVertex = 1, mesh % nVertices

       stressDivergenceU(iVertex) = 0.0_RKIND
       stressDivergenceV(iVertex) = 0.0_RKIND

       if (solveVelocity(iVertex) == 1) then

          ! loop over surrounding cells
          do iCellOnVertex = 1, mesh % vertexDegree
             
             ! get the cell number of this cell
             iCell = mesh % cellsOnVertex % array(iCellOnVertex, iVertex)

             nEdgesOnCell = mesh % nEdgesOnCell % array(iCell)

             alphaPWL = 1.0_RKIND / real(nEdgesOnCell, RKIND)
             
             ! get the vertexOnCell number of the iVertex velocity point from cell iCell
             jVertexOnCell = cellVerticesAtVertex(iCellOnVertex,iVertex)

             stressDivergenceUCell = 0.0_RKIND
             stressDivergenceVCell = 0.0_RKIND

             ! loop over the vertices of the surrounding cell
             do iVertexOnCell = 1, nEdgesOnCell

                ! loop over subcells
                do iSubCell = 1, nEdgesOnCell

                   stressDivergenceUCell = stressDivergenceUCell + &
                        stress11(iVertexOnCell,iCell) * basisSubArea(iSubCell,iCell) * &
                                                        pwl_basis_gradient(nEdgesOnCell, basisGradientU(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell) + &
                        stress12(iVertexOnCell,iCell) * basisSubArea(iSubCell,iCell) * &
                                                        pwl_basis_gradient(nEdgesOnCell, basisGradientV(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell)
                   
                   stressDivergenceVCell = stressDivergenceVCell + &
                        stress22(iVertexOnCell,iCell) * basisSubArea(iSubCell,iCell) * &
                                                        pwl_basis_gradient(nEdgesOnCell, basisGradientV(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell) + &
                        stress12(iVertexOnCell,iCell) * basisSubArea(iSubCell,iCell) * &
                                                        pwl_basis_gradient(nEdgesOnCell, basisGradientU(:,:,iCell), alphaPWL, iVertexOnCell, iSubCell)

                enddo ! iSubCell
                
             enddo ! jVertex

             stressDivergenceU(iVertex) = stressDivergenceU(iVertex) - stressDivergenceUCell
             stressDivergenceV(iVertex) = stressDivergenceV(iVertex) - stressDivergenceVCell
             
          enddo ! iCellOnVertex

          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) / mesh % areaTriangle % array(iVertex)
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) / mesh % areaTriangle % array(iVertex)

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
