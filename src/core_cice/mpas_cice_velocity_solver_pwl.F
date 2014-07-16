!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_velocity_solver_pwl
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_velocity_solver_pwl

  use mpas_grid_types

  implicit none

  private
  save

  public :: &
       cice_init_velocity_solver_pwl, &
       cice_internal_stress_pwl, &
       cice_strain_tensor_pwl, &
       cice_stress_divergence_pwl

contains

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_velocity_solver_pwl
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_init_velocity_solver_pwl(&
       mesh, &
       velocity_variational, &
       velocity_pwl, &
       boundary)!{{{

    use cice_mesh, only: &
         cice_cell_vertices_at_vertex

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    type(MPAS_pool_type), pointer :: &
         velocity_variational, &   !< Input/Output: 
         boundary, & !< Input/Output: 
         velocity_pwl      !< Input/Output: 

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex

    real(kind=RKIND), dimension(:,:), pointer :: &
         basisSubArea

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisGradientU, &
         basisGradientV

    call MPAS_pool_get_array(velocity_variational, "cellVerticesAtVertex", cellVerticesAtVertex)
    call MPAS_pool_get_array(velocity_pwl, "basisSubArea", basisSubArea)
    call MPAS_pool_get_array(velocity_pwl, "basisGradientU", basisGradientU)
    call MPAS_pool_get_array(velocity_pwl, "basisGradientV", basisGradientV)

    call cice_cell_vertices_at_vertex(&
         mesh, &
         cellVerticesAtVertex)

    call init_velocity_solver_pwl_basis(mesh, &
         basisSubArea, &
         basisGradientU, &
         basisGradientV)

  end subroutine cice_init_velocity_solver_pwl!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_velocity_solver_pwl_basis
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_velocity_solver_pwl_basis(&
       mesh, &
       basisSubArea, &
       basisGradientU, &
       basisGradientV)!{{{

    use cice_numerics, only: &
         cice_solve_linear_basis_system

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         basisSubArea !< Output: 
         
    real(kind=RKIND), dimension(:,:,:), intent(out) :: &    
         basisGradientU, & !< Output: 
         basisGradientV    !< Output: 

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
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

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

          call cice_solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

          basisGradientU(1,iEdgeOnCell,iCell) = solutionVector(1)
          basisGradientV(1,iEdgeOnCell,iCell) = solutionVector(2)

          ! second basis
          rightHandSide(1) = 0.0_RKIND
          rightHandSide(2) = 1.0_RKIND
          rightHandSide(3) = 0.0_RKIND

          call cice_solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

          basisGradientU(2,iEdgeOnCell,iCell) = solutionVector(1)
          basisGradientV(2,iEdgeOnCell,iCell) = solutionVector(2)

       enddo ! iEdgeOnCell

    enddo ! iCell

  end subroutine init_velocity_solver_pwl_basis!}}}

!-----------------------------------------------------------------------
! time step
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_internal_stress_pwl
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_internal_stress_pwl(&
       mesh, &
       velocity_variational, &
       velocity_pwl, &
       uVelocity, &
       vVelocity, &
       icePressure, &
       replacementPressure, &
       dtElastic, &
       evpDamping, &
       stressDivergenceU, &
       stressDivergenceV, &
       solveStress, &
       solveVelocity)!{{{

    use cice_velocity_solver_variational, only: &
         cice_stress_tensor_variational

    type (MPAS_pool_type), pointer :: &
         mesh, &   !< Input/Output: 
         velocity_variational, & !< Input/Output: 
         velocity_pwl    !< Input/Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, & !< Input: 
         vVelocity    !< Input:   

    real(kind=RKIND), dimension(:), intent(inout) :: &
         icePressure, &      !< Input/Output: 
         replacementPressure !< Input/Output: 

    real(kind=RKIND), intent(in) :: &
         dtElastic !< Input: 

    logical, intent(in) :: &
         evpDamping !< Input: 

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, & !< Output: 
         stressDivergenceV    !< Output: 

    integer, dimension(:), intent(in) :: &
         solveStress, & !< Input: 
         solveVelocity  !< Input: 

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

    call MPAS_pool_get_array(velocity_variational, "strain11", strain11)
    call MPAS_pool_get_array(velocity_variational, "strain22", strain22)
    call MPAS_pool_get_array(velocity_variational, "strain12", strain12)
    call MPAS_pool_get_array(velocity_variational, "stress11", stress11)
    call MPAS_pool_get_array(velocity_variational, "stress22", stress22)
    call MPAS_pool_get_array(velocity_variational, "stress12", stress12)
    call MPAS_pool_get_array(velocity_variational, "cellVerticesAtVertex", cellVerticesAtVertex)

    call MPAS_pool_get_array(velocity_pwl, "basisGradientU", basisGradientUPWL)
    call MPAS_pool_get_array(velocity_pwl, "basisGradientV", basisGradientVPWL)
    call MPAS_pool_get_array(velocity_pwl, "basisSubArea", basisSubArea)

    call cice_strain_tensor_pwl(mesh, &
         strain11,          &
         strain22,          &
         strain12,          &
         uVelocity,         &
         vVelocity,         &
         basisGradientUPWL, &
         basisGradientVPWL, &
         solveStress)
    
    call cice_stress_tensor_variational(mesh, &
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
    
    call cice_stress_divergence_pwl(mesh, &
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

  end subroutine cice_internal_stress_pwl!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_strain_tensor_pwl
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------  

  subroutine cice_strain_tensor_pwl(&
       mesh, &
       strain11, &
       strain22, &
       strain12, &
       uVelocity, &
       vVelocity, &
       basisGradientU, &
       basisGradientV, &
       solveStress)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         strain11, & !< Output: 
         strain22, & !< Output: 
         strain12    !< Output: 
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, & !< Input: 
         vVelocity    !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         basisGradientU, & !< Input: 
         basisGradientV    !< Input: 
    
    integer, dimension(:), intent(in) :: &
         solveStress !< Input: 

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
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

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

  end subroutine cice_strain_tensor_pwl!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_stress_divergence_pwl
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_stress_divergence_pwl(&
       mesh, &
       stressDivergenceU, &
       stressDivergenceV, &
       stress11, &
       stress22, &
       stress12, &
       basisGradientU, &
       basisGradientV, &
       basisSubArea, &
       cellVerticesAtVertex, &
       solveVelocity)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, & !< Output: 
         stressDivergenceV    !< Output: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &    
         stress11, & !< Input: 
         stress22, & !< Input: 
         stress12    !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &   
         basisGradientU, & !< Input: 
         basisGradientV    !< Input: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         basisSubArea !< Input: 

    integer, dimension(:,:), intent(in) :: &
         cellVerticesAtVertex !< Input: 

    integer, dimension(:), intent(in) :: &
         solveVelocity !< Input: 

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
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)

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
 
  end subroutine cice_stress_divergence_pwl!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  pwl_basis_gradient
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function pwl_basis_gradient(&
       nEdgesOnCell, &
       basisGradient, &
       alphaPWL, &
       basisVertexOnCell, &
       iSubCell) &
       result(grad)!{{{
    
    integer, intent(in) :: &
         nEdgesOnCell, &      !< Input: 
         basisVertexOnCell, & !< Input: basis function vertex
         iSubCell             !< Input: subcell to calculate

    real(kind=RKIND), intent(in) :: &
         alphaPWL !< Input: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         basisGradient !< Input: 

    real(kind=RKIND) :: grad
    
    grad = basisGradient(1,iSubCell) * &
           (alphaPWL + merge(1.0_RKIND, 0.0_RKIND, iSubCell == basisVertexOnCell)) + &
           basisGradient(2,iSubCell) * &
           (alphaPWL + merge(1.0_RKIND, 0.0_RKIND, iSubCell == wrapped_index(basisVertexOnCell - 1, nEdgesOnCell)))

  end function pwl_basis_gradient!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wrapped_index
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function wrapped_index(&
       input, &
       nelements) &
       result(output)!{{{

    integer, intent(in) :: &
         input, &  !< Input: 
         nelements !< Input: 

    integer :: output

    output = mod(input - 1, nelements) + 1

  end function wrapped_index!}}}

!-----------------------------------------------------------------------

end module cice_velocity_solver_pwl
