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
       cice_init_velocity_solver_pwl

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
       boundary, &
       rotateCartesianGrid, &
       includeMetricTerms)!{{{

    use cice_mesh, only: &
         cice_cell_vertices_at_vertex

    use cice_velocity_solver_variational_shared, only: &
         cice_calc_local_coords, &
         cice_calc_variational_metric_terms
    
    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    type(MPAS_pool_type), pointer :: &
         velocity_variational, & !< Input/Output: 
         boundary                !< Input/Output: 

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input: 
         includeMetricTerms     !< Input: 

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex

    integer, pointer :: &
         nCells, &
         maxEdges

    real(kind=RKIND), dimension(:), pointer :: &
         tanLatVertexRotatedOverRadius

    real(kind=RKIND), dimension(:,:), allocatable :: &
         xLocal, &
         yLocal

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisGradientU, &
         basisGradientV, &
         basisIntegralsMetric, &
         basisIntegralsU, &
         basisIntegralsV

    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "maxEdges", maxEdges)

    call MPAS_pool_get_array(velocity_variational, "cellVerticesAtVertex", cellVerticesAtVertex)
    call MPAS_pool_get_array(velocity_variational, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)
    call MPAS_pool_get_array(velocity_variational, "basisGradientU", basisGradientU)
    call MPAS_pool_get_array(velocity_variational, "basisGradientV", basisGradientV)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsMetric", basisIntegralsMetric)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsU", basisIntegralsU)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsV", basisIntegralsV)

    allocate(xLocal(maxEdges,nCells))
    allocate(yLocal(maxEdges,nCells))

    call cice_calc_local_coords(&
         mesh, &
         xLocal, &
         yLocal, &
         rotateCartesianGrid)

    call cice_calc_variational_metric_terms(&
         mesh, &
         tanLatVertexRotatedOverRadius, &
         rotateCartesianGrid, &
         includeMetricTerms)

    call cice_cell_vertices_at_vertex(&
         mesh, &
         cellVerticesAtVertex)

    call init_velocity_solver_pwl_basis(&
         mesh, &
         basisGradientU, &
         basisGradientV, &
         basisIntegralsMetric, &
         basisIntegralsU, &
         basisIntegralsV, &
         xLocal, &
         yLocal)

    deallocate(xLocal)
    deallocate(yLocal)

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
       basisGradientU, &
       basisGradientV, &
       basisIntegralsMetric, &
       basisIntegralsU, &
       basisIntegralsV, &
       xLocal, &
       yLocal)!{{{

    use cice_numerics, only: &
         cice_solve_linear_basis_system

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &    
         basisGradientU, &       !< Output: 
         basisGradientV, &       !< Output: 
         basisIntegralsMetric, & !< Output:
         basisIntegralsU, &      !< Output:
         basisIntegralsV         !< Output:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         xLocal, & !< Input: 
         yLocal    !< Input: 

    real(kind=RKIND) :: &
         xPWLCentre, &
         yPWLCentre, &
         alphaPWL, &
         a, b, c, s, &
         basisIntegral, &
         basisIntegralsMetricSubCell

    integer :: &
         iCell, &
         iVertexOnCell, &
         iEdgeOnCell, &
         iEdge, &
         iVertexOnCell1, &
         iVertexOnCell2, &
         iSubCell, &
         jEdgeOnCell, &
         jEdgeOnCell2

    real(kind=RKIND), dimension(3,3) :: &
         leftMatrix

    real(kind=RKIND), dimension(3) :: &
         rightHandSide, &
         solutionVector

    integer, pointer :: &
         nCells, &
         maxEdges

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell, &
         edgesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         dvEdge

    real(kind=RKIND), dimension(:,:), allocatable :: &
         subBasisGradientU, &
         subBasisGradientV, &
         subCellgradientU, &
         subCellgradientV

    real(kind=RKIND), dimension(:), allocatable :: &
         basisSubArea

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "maxEdges", maxEdges)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "dvEdge", dvEdge)

    allocate(subBasisGradientU(maxEdges,3))
    allocate(subBasisGradientV(maxEdges,3))
    allocate(subCellgradientU(maxEdges,maxEdges))
    allocate(subCellgradientV(maxEdges,maxEdges))
    allocate(basisSubArea(maxEdges))

    ! loop over cells
    do iCell = 1, nCells
       
       alphaPWL = 1.0_RKIND / real(nEdgesOnCell(iCell),RKIND)

       ! determine cell centre for piecewise linear basis
       xPWLCentre = 0.0_RKIND
       yPWLCentre = 0.0_RKIND

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          xPWLCentre = xPWLCentre + alphaPWL * xLocal(iVertexOnCell,iCell)
          yPWLCentre = yPWLCentre + alphaPWL * yLocal(iVertexOnCell,iCell)

       enddo ! iVertexOnCell

       ! calculate the area of the subcells
       do iSubCell = 1, nEdgesOnCell(iCell)

          iEdge = edgesOnCell(iSubCell,iCell)
          iVertexOnCell1 = iSubCell
          iVertexOnCell2 = wrapped_index(iSubCell + 1, nEdgesOnCell(iCell))

          c = dvEdge(iEdge)
          a = sqrt((xLocal(iVertexOnCell1,iCell) - xPWLCentre)**2 + &
                   (yLocal(iVertexOnCell1,iCell) - yPWLCentre)**2)
          b = sqrt((xLocal(iVertexOnCell2,iCell) - xPWLCentre)**2 + &
                   (yLocal(iVertexOnCell2,iCell) - yPWLCentre)**2)

          s = (a + b + c) * 0.5_RKIND

          ! Heron's formula 
          basisSubArea(iSubCell) = sqrt(s * (s-a) * (s-b) * (s-c))

       enddo ! iSubCell

       ! calculate the linear basis on the sub triangle
       do iSubCell = 1, nEdgesOnCell(iCell)

          iVertexOnCell1 = iSubCell
          iVertexOnCell2 = wrapped_index(iSubCell + 1, nEdgesOnCell(iCell))

          ! set up left hand matrix
          leftMatrix(1,1) = xLocal(iVertexOnCell1,iCell) - xPWLCentre
          leftMatrix(1,2) = yLocal(iVertexOnCell1,iCell) - yPWLCentre
          leftMatrix(1,3) = 1.0_RKIND
          
          leftMatrix(2,1) = xLocal(iVertexOnCell2,iCell) - xPWLCentre
          leftMatrix(2,2) = yLocal(iVertexOnCell2,iCell) - yPWLCentre
          leftMatrix(2,3) = 1.0_RKIND
          
          leftMatrix(3,1) = 0.0_RKIND
          leftMatrix(3,2) = 0.0_RKIND
          leftMatrix(3,3) = 1.0_RKIND

          ! first basis
          rightHandSide(1) = 1.0_RKIND
          rightHandSide(2) = 0.0_RKIND
          rightHandSide(3) = 0.0_RKIND

          call cice_solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

          subBasisGradientU(1,iSubCell) = solutionVector(1)
          subBasisGradientV(1,iSubCell) = solutionVector(2)

          ! second basis
          rightHandSide(1) = 0.0_RKIND
          rightHandSide(2) = 1.0_RKIND
          rightHandSide(3) = 0.0_RKIND

          call cice_solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

          subBasisGradientU(2,iSubCell) = solutionVector(1)
          subBasisGradientV(2,iSubCell) = solutionVector(2)

          ! third basis
          subBasisGradientU(3,iSubCell) = 1.0_RKIND - subBasisGradientU(1,iSubCell) - subBasisGradientU(2,iSubCell)
          subBasisGradientV(3,iSubCell) = 1.0_RKIND - subBasisGradientV(1,iSubCell) - subBasisGradientV(2,iSubCell)

       enddo ! iSubCell

       ! use the linear sub area basis to calculate the PWL basis
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          ! loop over subcells
          do iSubCell = 1, nEdgesOnCell(iCell)

             ! array (index of the basis vertex, subarea value)
             subCellGradientU(iEdgeOnCell,iSubCell) = subBasisGradientU(3,iSubCell) * alphaPWL
             subCellGradientV(iEdgeOnCell,iSubCell) = subBasisGradientV(3,iSubCell) * alphaPWL

             if (iSubCell == iEdgeOnCell) then

                subCellGradientU(iEdgeOnCell,iSubCell) = subCellGradientU(iEdgeOnCell,iSubCell) + subBasisGradientU(1,iSubCell)
                subCellGradientV(iEdgeOnCell,iSubCell) = subCellGradientV(iEdgeOnCell,iSubCell) + subBasisGradientV(1,iSubCell)

             else if (iSubCell == wrapped_index(iEdgeOnCell - 1, nEdgesOnCell(iCell))) then

                subCellGradientU(iEdgeOnCell,iSubCell) = subCellGradientU(iEdgeOnCell,iSubCell) + subBasisGradientU(2,iSubCell)
                subCellGradientV(iEdgeOnCell,iSubCell) = subCellGradientV(iEdgeOnCell,iSubCell) + subBasisGradientV(2,iSubCell)

             endif

          enddo ! iSubCell
          
       enddo ! iEdgeOnCell

       ! calculate the gradients at the cell corners
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          do jEdgeOnCell = 1, nEdgesOnCell(iCell)

             jEdgeOnCell2 = wrapped_index(jEdgeOnCell - 1, nEdgesOnCell(iCell))

             basisGradientU(iEdgeOnCell,jEdgeOnCell,iCell) = 0.5_RKIND * &
                  (subCellGradientU(iEdgeOnCell,jEdgeOnCell) + subCellGradientU(iEdgeOnCell,jEdgeOnCell2))
             basisGradientV(iEdgeOnCell,jEdgeOnCell,iCell) = 0.5_RKIND * &
                  (subCellGradientV(iEdgeOnCell,jEdgeOnCell) + subCellGradientV(iEdgeOnCell,jEdgeOnCell2))

          enddo ! jEdgeOnCell

       enddo ! iEdgeOnCell

       ! calculate the basis integrals
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)
          
          do jEdgeOnCell = 1, nEdgesOnCell(iCell)

             basisIntegralsU(iEdgeOnCell,jEdgeOnCell,iCell) = 0.0_RKIND
             basisIntegralsV(iEdgeOnCell,jEdgeOnCell,iCell) = 0.0_RKIND

             do iSubCell = 1, nEdgesOnCell(iCell)

                if (iSubCell == iEdgeOnCell .or. iSubCell == wrapped_index(iEdgeOnCell - 1, nEdgesOnCell(iCell))) then
                   basisIntegral = ((alphaPWL + 1) * basisSubArea(iSubCell)) / 3.0_RKIND
                else
                   basisIntegral = ( alphaPWL      * basisSubArea(iSubCell)) / 3.0_RKIND
                endif

                basisIntegralsU(iEdgeOnCell,jEdgeOnCell,iCell) = basisIntegralsU(iEdgeOnCell,jEdgeOnCell,iCell) + &
                     subCellGradientU(jEdgeOnCell,iSubCell) * basisIntegral

                basisIntegralsV(iEdgeOnCell,jEdgeOnCell,iCell) = basisIntegralsV(iEdgeOnCell,jEdgeOnCell,iCell) + &
                     subCellGradientV(jEdgeOnCell,iSubCell) * basisIntegral

             enddo ! iSubCell

          enddo ! jEdgeOnCell

       enddo ! iEdgeOnCell

       ! basis integrals for the metric terms
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)
          do jEdgeOnCell = 1, nEdgesOnCell(iCell)

             basisIntegralsMetric(iEdgeOnCell,jEdgeOnCell,iCell) = 0.0_RKIND

             do iSubCell = 1, nEdgesOnCell(iCell)

                basisIntegralsMetricSubCell = alphaPWL**2

                if ( (iSubCell == iEdgeOnCell .and. &
                     iSubCell == jEdgeOnCell) .or. &
                     (iSubCell == wrapped_index(iEdgeOnCell - 1, nEdgesOnCell(iCell)) .and. &
                      iSubCell == wrapped_index(jEdgeOnCell - 1, nEdgesOnCell(iCell))) ) then

                   basisIntegralsMetricSubCell = 2.0_RKIND * (alphaPWL + 1)**2

                endif

                if ( (iSubCell == iEdgeOnCell .and. &
                      iSubCell == wrapped_index(jEdgeOnCell - 1, nEdgesOnCell(iCell))) .or. &
                     (iSubCell == wrapped_index(iEdgeOnCell - 1, nEdgesOnCell(iCell)) .and. &
                      iSubCell == jEdgeOnCell) ) then

                   basisIntegralsMetricSubCell = 2.0_RKIND * alphaPWL**2 + 2.0_RKIND * alphaPWL + 1.0_RKIND

                endif

                basisIntegralsMetricSubCell = basisIntegralsMetricSubCell * &
                     basisSubArea(iSubCell) / 12.0_RKIND

                basisIntegralsMetric(iEdgeOnCell,jEdgeOnCell,iCell) = basisIntegralsMetric(iEdgeOnCell,jEdgeOnCell,iCell) + &
                     basisIntegralsMetricSubCell

             enddo ! iSubCell

          enddo ! jEdgeOnCell
       enddo ! iSubCell

    enddo ! iCell

    deallocate(subBasisGradientU)
    deallocate(subBasisGradientV)
    deallocate(subCellgradientU)
    deallocate(subCellgradientV)
    deallocate(basisSubArea)

  end subroutine init_velocity_solver_pwl_basis!}}}

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
