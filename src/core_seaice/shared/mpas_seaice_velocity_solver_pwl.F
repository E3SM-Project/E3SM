!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_velocity_solver_pwl
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_velocity_solver_pwl

  use mpas_derived_types
  use mpas_pool_routines

  implicit none

  private
  save

  public :: &
       seaice_init_velocity_solver_pwl

contains

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_velocity_solver_pwl
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_velocity_solver_pwl(&
       mesh, &
       velocity_variational, &
       boundary, &
       rotateCartesianGrid, &
       includeMetricTerms)!{{{

    use seaice_mesh, only: &
         seaice_cell_vertices_at_vertex

    use seaice_velocity_solver_variational_shared, only: &
         seaice_calc_local_coords, &
         seaice_calc_variational_metric_terms

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

    integer :: iCell, i1, i2

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "maxEdges", maxEdges)

    call MPAS_pool_get_array(velocity_variational, "cellVerticesAtVertex", cellVerticesAtVertex)
    call MPAS_pool_get_array(velocity_variational, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)
    call MPAS_pool_get_array(velocity_variational, "basisGradientU", basisGradientU)
    call MPAS_pool_get_array(velocity_variational, "basisGradientV", basisGradientV)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsMetric", basisIntegralsMetric)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsU", basisIntegralsU)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsV", basisIntegralsV)

    allocate(xLocal(maxEdges,nCells))
    allocate(yLocal(maxEdges,nCells))

    call seaice_calc_local_coords(&
         mesh, &
         xLocal, &
         yLocal, &
         rotateCartesianGrid)

    call seaice_calc_variational_metric_terms(&
         mesh, &
         tanLatVertexRotatedOverRadius, &
         rotateCartesianGrid, &
         includeMetricTerms)

    call seaice_cell_vertices_at_vertex(&
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

  end subroutine seaice_init_velocity_solver_pwl!}}}

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

    use seaice_numerics, only: &
         seaice_solve_linear_basis_system

    use seaice_velocity_solver_variational_shared, only: &
         seaice_wrapped_index

    type(MPAS_pool_type), pointer :: &
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
         basisIntegralsMetricSubCell, &
         basisSubAreaSum

    integer :: &
         iCell, &
         iVertexOnCell, &
         iEdgeOnCell, &
         iEdge, &
         iVertexOnCell1, &
         iVertexOnCell2, &
         iSubCell, &
         iBasisVertex, &
         iGradientVertex, &
         iSubCell1, &
         iSubCell2, &
         iStressVertex, &
         iVelocityVertex, &
         subCellTypeStress, &
         subCellTypeVelocity

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
         dvEdge, &
         areaCell

    real(kind=RKIND), dimension(:,:), allocatable :: &
         subBasisGradientU, &
         subBasisGradientV, &
         subBasisConstant, &
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
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    allocate(subBasisGradientU(maxEdges,3))
    allocate(subBasisGradientV(maxEdges,3))
    allocate(subBasisConstant(maxEdges,3))
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
       basisSubAreaSum = 0.0_RKIND

       do iSubCell = 1, nEdgesOnCell(iCell)

          iEdge = edgesOnCell(iSubCell,iCell)
          iVertexOnCell1 = iSubCell
          iVertexOnCell2 = seaice_wrapped_index(iSubCell + 1, nEdgesOnCell(iCell))

          c = dvEdge(iEdge)
          a = sqrt((xLocal(iVertexOnCell1,iCell) - xPWLCentre)**2 + &
                   (yLocal(iVertexOnCell1,iCell) - yPWLCentre)**2)
          b = sqrt((xLocal(iVertexOnCell2,iCell) - xPWLCentre)**2 + &
                   (yLocal(iVertexOnCell2,iCell) - yPWLCentre)**2)

          s = (a + b + c) * 0.5_RKIND

          ! Heron's formula
          basisSubArea(iSubCell) = sqrt(s * (s-a) * (s-b) * (s-c))

          basisSubAreaSum = basisSubAreaSum + basisSubArea(iSubCell)

       enddo ! iSubCell

       ! ensure sum of subareas equals the area of the cell
       basisSubArea(:) = basisSubArea(:) * (areaCell(iCell) / basisSubAreaSum)

       ! calculate the linear basis on the sub triangle
       do iSubCell = 1, nEdgesOnCell(iCell)

          iVertexOnCell1 = iSubCell
          iVertexOnCell2 = seaice_wrapped_index(iSubCell + 1, nEdgesOnCell(iCell))

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

          call seaice_solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

          subBasisGradientU(iSubCell,1) = solutionVector(1)
          subBasisGradientV(iSubCell,1) = solutionVector(2)
          subBasisConstant(iSubCell,1)  = solutionVector(3)

          ! second basis
          rightHandSide(1) = 0.0_RKIND
          rightHandSide(2) = 1.0_RKIND
          rightHandSide(3) = 0.0_RKIND

          call seaice_solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

          subBasisGradientU(iSubCell,2) = solutionVector(1)
          subBasisGradientV(iSubCell,2) = solutionVector(2)
          subBasisConstant(iSubCell,2)  = solutionVector(3)

          ! third basis
          subBasisGradientU(iSubCell,3) = -subBasisGradientU(iSubCell,1) - subBasisGradientU(iSubCell,2)
          subBasisGradientV(iSubCell,3) = -subBasisGradientV(iSubCell,1) - subBasisGradientV(iSubCell,2)
          subBasisConstant(iSubCell,3)  = 1.0_RKIND - subBasisConstant(iSubCell,1) - subBasisConstant(iSubCell,2)

       enddo ! iSubCell

       ! use the linear sub area basis to calculate the PWL basis
       do iBasisVertex = 1, nEdgesOnCell(iCell)

          ! loop over subcells
          do iSubCell = 1, nEdgesOnCell(iCell)

             ! array (index of the basis vertex, subarea value)
             subCellGradientU(iBasisVertex,iSubCell) = subBasisGradientU(iSubCell,3) * alphaPWL
             subCellGradientV(iBasisVertex,iSubCell) = subBasisGradientV(iSubCell,3) * alphaPWL

             if (iSubCell == iBasisVertex) then

                subCellGradientU(iBasisVertex,iSubCell) = subCellGradientU(iBasisVertex,iSubCell) + subBasisGradientU(iSubCell,1)
                subCellGradientV(iBasisVertex,iSubCell) = subCellGradientV(iBasisVertex,iSubCell) + subBasisGradientV(iSubCell,1)

             else if (iSubCell == seaice_wrapped_index(iBasisVertex - 1, nEdgesOnCell(iCell))) then

                subCellGradientU(iBasisVertex,iSubCell) = subCellGradientU(iBasisVertex,iSubCell) + subBasisGradientU(iSubCell,2)
                subCellGradientV(iBasisVertex,iSubCell) = subCellGradientV(iBasisVertex,iSubCell) + subBasisGradientV(iSubCell,2)

             endif

          enddo ! iSubCell

       enddo ! iEdgeOnCell

       ! calculate the gradients at the cell corners
       do iBasisVertex = 1, nEdgesOnCell(iCell)

          do iGradientVertex = 1, nEdgesOnCell(iCell)

             iSubCell1 = iGradientVertex
             iSubCell2 = seaice_wrapped_index(iGradientVertex - 1, nEdgesOnCell(iCell))

             basisGradientU(iBasisVertex,iGradientVertex,iCell) = &
                  0.5_RKIND * (subCellGradientU(iBasisVertex,iSubCell1) + subCellGradientU(iBasisVertex,iSubCell2))
             basisGradientV(iBasisVertex,iGradientVertex,iCell) = &
                  0.5_RKIND * (subCellGradientV(iBasisVertex,iSubCell1) + subCellGradientV(iBasisVertex,iSubCell2))

          enddo ! iGradientVertex

       enddo ! iBasisVertex

       ! calculate the basis integrals
       do iStressVertex = 1, nEdgesOnCell(iCell)
          do iVelocityVertex = 1, nEdgesOnCell(iCell)

             basisIntegralsU(iStressVertex,iVelocityVertex,iCell) = 0.0_RKIND
             basisIntegralsV(iStressVertex,iVelocityVertex,iCell) = 0.0_RKIND

             do iSubCell = 1, nEdgesOnCell(iCell)

                if (iSubCell == iStressVertex .or. iSubCell == seaice_wrapped_index(iStressVertex - 1, nEdgesOnCell(iCell))) then
                   basisIntegral = ((alphaPWL + 1) * basisSubArea(iSubCell)) / 3.0_RKIND
                else
                   basisIntegral = ( alphaPWL      * basisSubArea(iSubCell)) / 3.0_RKIND
                endif

                basisIntegralsU(iStressVertex,iVelocityVertex,iCell) = basisIntegralsU(iStressVertex,iVelocityVertex,iCell) + &
                     subCellGradientU(iVelocityVertex,iSubCell) * basisIntegral

                basisIntegralsV(iStressVertex,iVelocityVertex,iCell) = basisIntegralsV(iStressVertex,iVelocityVertex,iCell) + &
                     subCellGradientV(iVelocityVertex,iSubCell) * basisIntegral

             enddo ! iSubCell

          enddo ! iVelocityVertex
       enddo ! iStressVertex

       ! basis integrals for the metric terms
       do iStressVertex = 1, nEdgesOnCell(iCell)
          do iVelocityVertex = 1, nEdgesOnCell(iCell)

             basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) = 0.0_RKIND

             do iSubCell = 1, nEdgesOnCell(iCell)

                ! determine stress subcell type
                if (iSubCell == iStressVertex) then
                   subCellTypeStress = 1
                else if (iSubCell == seaice_wrapped_index(iStressVertex - 1, nEdgesOnCell(iCell))) then
                   subCellTypeStress = 2
                else
                   subCellTypeStress = 3
                endif

                ! determine velocity subcell type
                if (iSubCell == iVelocityVertex) then
                   subCellTypeVelocity = 1
                else if (iSubCell == seaice_wrapped_index(iVelocityVertex - 1, nEdgesOnCell(iCell))) then
                   subCellTypeVelocity = 2
                else
                   subCellTypeVelocity = 3
                endif

                ! set the subcell integral value
                if      ((subCellTypeStress == 1 .and. subCellTypeVelocity == 1) .or. &
                         (subCellTypeStress == 2 .and. subCellTypeVelocity == 2)) then

                   basisIntegralsMetricSubCell = 2.0_RKIND * alphaPWL**2 + 2.0_RKIND * alphaPWL + 2.0_RKIND

                else if ((subCellTypeStress == 1 .and. subCellTypeVelocity == 2) .or. &
                         (subCellTypeStress == 2 .and. subCellTypeVelocity == 1)) then

                   basisIntegralsMetricSubCell = 2.0_RKIND * alphaPWL**2 + 2.0_RKIND * alphaPWL + 1.0_RKIND

                else if ((subCellTypeStress == 1 .and. subCellTypeVelocity == 3) .or. &
                         (subCellTypeStress == 3 .and. subCellTypeVelocity == 1) .or. &
                         (subCellTypeStress == 2 .and. subCellTypeVelocity == 3) .or. &
                         (subCellTypeStress == 3 .and. subCellTypeVelocity == 2)) then

                   basisIntegralsMetricSubCell = 2.0_RKIND * alphaPWL**2 + alphaPWL

                else if  (subCellTypeStress == 3 .and. subCellTypeVelocity == 3) then

                   basisIntegralsMetricSubCell = 2.0_RKIND * alphaPWL**2

                end if

                basisIntegralsMetricSubCell = basisIntegralsMetricSubCell * &
                     basisSubArea(iSubCell) / 12.0_RKIND

                basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) = &
                     basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) + &
                     basisIntegralsMetricSubCell

             enddo ! iSubCell

          enddo ! iVelocityVertex
       enddo ! iStressVertex

    enddo ! iCell

    deallocate(subBasisGradientU)
    deallocate(subBasisGradientV)
    deallocate(subBasisConstant)
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

    use seaice_velocity_solver_variational_shared, only: &
         seaice_wrapped_index

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
           (alphaPWL + merge(1.0_RKIND, 0.0_RKIND, iSubCell == seaice_wrapped_index(basisVertexOnCell - 1, nEdgesOnCell)))

  end function pwl_basis_gradient!}}}

!-----------------------------------------------------------------------

end module seaice_velocity_solver_pwl
