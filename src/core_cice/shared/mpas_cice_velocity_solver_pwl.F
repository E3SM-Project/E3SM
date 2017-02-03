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

  use mpas_derived_types
  use mpas_pool_routines

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

    use cice_debug, only: &
         cice_plot_cells_planar

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

    use cice_velocity_solver_variational_shared, only: &
         cice_wrapped_index

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
          iVertexOnCell2 = cice_wrapped_index(iSubCell + 1, nEdgesOnCell(iCell))

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
          iVertexOnCell2 = cice_wrapped_index(iSubCell + 1, nEdgesOnCell(iCell))

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

          subBasisGradientU(iSubCell,1) = solutionVector(1)
          subBasisGradientV(iSubCell,1) = solutionVector(2)
          subBasisConstant(iSubCell,1)  = solutionVector(3)

          ! second basis
          rightHandSide(1) = 0.0_RKIND
          rightHandSide(2) = 1.0_RKIND
          rightHandSide(3) = 0.0_RKIND

          call cice_solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

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

             else if (iSubCell == cice_wrapped_index(iBasisVertex - 1, nEdgesOnCell(iCell))) then

                subCellGradientU(iBasisVertex,iSubCell) = subCellGradientU(iBasisVertex,iSubCell) + subBasisGradientU(iSubCell,2)
                subCellGradientV(iBasisVertex,iSubCell) = subCellGradientV(iBasisVertex,iSubCell) + subBasisGradientV(iSubCell,2)

             endif

          enddo ! iSubCell

       enddo ! iEdgeOnCell

       ! calculate the gradients at the cell corners
       do iBasisVertex = 1, nEdgesOnCell(iCell)

          do iGradientVertex = 1, nEdgesOnCell(iCell)

             iSubCell1 = iGradientVertex
             iSubCell2 = cice_wrapped_index(iGradientVertex - 1, nEdgesOnCell(iCell))

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

                if (iSubCell == iStressVertex .or. iSubCell == cice_wrapped_index(iStressVertex - 1, nEdgesOnCell(iCell))) then
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
                else if (iSubCell == cice_wrapped_index(iStressVertex - 1, nEdgesOnCell(iCell))) then
                   subCellTypeStress = 2
                else
                   subCellTypeStress = 3
                endif

                ! determine velocity subcell type
                if (iSubCell == iVelocityVertex) then
                   subCellTypeVelocity = 1
                else if (iSubCell == cice_wrapped_index(iVelocityVertex - 1, nEdgesOnCell(iCell))) then
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

    use cice_velocity_solver_variational_shared, only: &
         cice_wrapped_index

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
           (alphaPWL + merge(1.0_RKIND, 0.0_RKIND, iSubCell == cice_wrapped_index(basisVertexOnCell - 1, nEdgesOnCell)))

  end function pwl_basis_gradient!}}}

!-----------------------------------------------------------------------
! Hex specific plotting stuff
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_pwl
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 28 October 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine plot_pwl(&
       mesh, &
       iCell, &
       nEdgesOnCell, &
       iVertex, &
       subBasisGradientU, &
       subBasisGradientV, &
       subBasisConstant, &
       xLocal, &
       yLocal)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    integer, intent(in) :: &
         iCell, &        !< Input:
         nEdgesOnCell, & !< Input:
         iVertex         !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         subBasisGradientU, &
         subBasisGradientV, &
         subBasisConstant

    real(kind=RKIND), dimension(:), intent(in) :: &
         xLocal, & !< Input:
         yLocal    !< Input:

    real(kind=RKIND) :: &
         xmin, xmax, ymin, ymax

    real(kind=RKIND), dimension(2,2) :: &
         mapping

    integer :: &
         iSubTriangle, &
         i1, &
         i2, &
         iVertexOnCell, &
         iVertex2

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         xCell, &
         yCell

    ! init variables
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)

    xmin =  1e30
    xmax = -1e30
    ymin =  1e30
    ymax = -1e30

    do iVertexOnCell = 1, nEdgesOnCell

       iVertex2 = verticesOnCell(iVertexOnCell,iCell)

       xmin = min(xmin, xVertex(iVertex2)) - 1000.0_RKIND
       xmax = max(xmax, xVertex(iVertex2)) + 1000.0_RKIND
       ymin = min(ymin, yVertex(iVertex2)) - 1000.0_RKIND
       ymax = max(ymax, yVertex(iVertex2)) + 1000.0_RKIND

       open(55,file='test.txt')
       if (iVertexOnCell == iVertex) then
          write(55,*) xVertex(iVertex2), yVertex(iVertex2)
       endif
       close(55)

    enddo ! iVertex

    xmin = xCell(iCell) - 11000.0_RKIND
    xmax = xCell(iCell) + 11000.0_RKIND
    ymin = yCell(iCell) - 11000.0_RKIND
    ymax = yCell(iCell) + 11000.0_RKIND

    open(55,file="pwl.txt")

    write(55,fmt='(a,f10.2,a,f10.2,a)') "set xrange [",xmin,":",xmax,"]"
    write(55,fmt='(a,f10.2,a,f10.2,a)') "set yrange [",ymin,":",ymax,"]"

    ! loop over subtriangles
    do iSubTriangle = 1, nEdgesOnCell

       i1 = iSubTriangle
       i2 = iSubTriangle + 1
       if (i2 > nEdgesOnCell) i2 = i2 - nEdgesOnCell

       ! get the triangle mapping
       call get_triangle_mapping(mapping, &
                                 1.0_RKIND,0.0_RKIND,&
                                 0.0_RKIND,1.0_RKIND,&
                                 xLocal(i1),yLocal(i1),&
                                 xLocal(i2),yLocal(i2))

       ! plot from subtriangle
       call plot_subtriangle(&
            nEdgesOnCell, &
            iSubTriangle, &
            iVertex, &
            subBasisGradientU, &
            subBasisGradientV, &
            subBasisConstant, &
            mapping, &
            xCell(iCell), &
            yCell(iCell))

    enddo ! iSubTriangle

    close(55)

  end subroutine plot_pwl!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_subtriangle
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 28 October 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine plot_subtriangle(&
       nEdgesOnCell, &
       iSubTriangle, &
       iVertex, &
       subBasisGradientU, &
       subBasisGradientV, &
       subBasisConstant, &
       mapping, &
       x0, &
       y0)!{{{

    use cice_velocity_solver_variational_shared, only: &
         cice_wrapped_index

    integer, intent(in) :: &
         nEdgesOnCell, & !< Input:
         iSubTriangle, & !< Input:
         iVertex         !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         subBasisGradientU, & !< Input:
         subBasisGradientV, & !< Input:
         subBasisConstant     !< Input:

    real(kind=RKIND), dimension(2,2), intent(in) :: &
         mapping !< Input:

    real(kind=RKIND), intent(in) :: &
         x0, & !< Input:
         y0    !< Input:

    real(kind=RKIND) :: &
         u, v, &
         x, y, &
         pwl, &
         alpha

    integer :: &
         i, j, iObject

    real(kind=RKIND) :: &
         x1, x2, x3, x4, x5, &
         y1, y2, y3, y4, y5, &
         d, &
         dx, dy

    integer, parameter :: n = 50

    logical, parameter :: &
         lfinitedifference = .false.

    iObject = ((iSubTriangle - 1) * (n+1) * (n+2)) / 2

    d = 1.0_RKIND / real(n, RKIND)

    do i = 0, n
       do j = 0, n

          if (i<=n-j-1) then

             iObject = iObject + 1

             u = real(i,RKIND) / real(n,RKIND)
             v = real(j,RKIND) / real(n,RKIND)

             call use_triangle_mapping(x, y, u, v, mapping)

             alpha = 1.0_RKIND / real(nEdgesOnCell,RKIND)

             pwl = alpha * (subBasisGradientU(iSubTriangle,3) * x + &
                            subBasisGradientV(iSubTriangle,3) * y + &
                            subBasisConstant(iSubTriangle,3))

             write(*,*) iVertex, nEdgesOnCell, cice_wrapped_index(iVertex - 1, nEdgesOnCell)

             if (iSubTriangle == iVertex) then

                pwl = pwl + &
                     (subBasisGradientU(iSubTriangle,1) * x + &
                      subBasisGradientV(iSubTriangle,1) * y + &
                      subBasisConstant(iSubTriangle,1))

             else if (iSubTriangle == cice_wrapped_index(iVertex - 1, nEdgesOnCell)) then

                pwl = pwl + &
                     (subBasisGradientU(iSubTriangle,2) * x + &
                      subBasisGradientV(iSubTriangle,2) * y + &
                      subBasisConstant(iSubTriangle,2))

             endif

             if (i<n-j-1) then

                call use_triangle_mapping(x1, y1, u,   v,   mapping)
                call use_triangle_mapping(x2, y2, u,   v+d, mapping)
                call use_triangle_mapping(x3, y3, u+d, v+d, mapping)
                call use_triangle_mapping(x4, y4, u+d, v,   mapping)
                x5 = x1
                y5 = y1

                x1 = x1 + x0 ; y1 = y1 + y0
                x2 = x2 + x0 ; y2 = y2 + y0
                x3 = x3 + x0 ; y3 = y3 + y0
                x4 = x4 + x0 ; y4 = y4 + y0
                x5 = x5 + x0 ; y5 = y5 + y0

                write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') &
                     "set object ",iObject," polygon from ",&
                     x1,",",y1," to ",x2,",",y2," to ",x3,",",y3," to ",x4,",",y4," to ",x5,",",y5

                write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', pwl, ' fillstyle solid'

             else

                call use_triangle_mapping(x1, y1, u,   v,   mapping)
                call use_triangle_mapping(x2, y2, u,   v+d, mapping)
                call use_triangle_mapping(x4, y4, u+d, v,   mapping)
                x5 = x1
                y5 = y1

                x1 = x1 + x0 ; y1 = y1 + y0
                x2 = x2 + x0 ; y2 = y2 + y0
                x4 = x4 + x0 ; y4 = y4 + y0
                x5 = x5 + x0 ; y5 = y5 + y0

                write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') &
                     "set object ",iObject," polygon from ",&
                     x1,",",y1," to ",x2,",",y2," to ",x4,",",y4," to ",x5,",",y5

                write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', pwl, ' fillstyle solid'

             endif

          endif

       enddo ! j
    enddo ! i

  end subroutine plot_subtriangle!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_triangle_mapping
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_triangle_mapping(&
       mapping, &
       x1, y1, &
       x2, y2, &
       u1, v1, &
       u2, v2)!{{{

    real(kind=RKIND), dimension(2,2), intent(out) :: &
         mapping !< Output:

    real(kind=RKIND), intent(in) :: &
         x1, & !< Input:
         y1, & !< Input:
         x2, & !< Input:
         y2, & !< Input:
         u1, & !< Input:
         v1, & !< Input:
         u2, & !< Input:
         v2    !< Input:

    mapping(1,1) = (u2*y1 - u1*y2) / (x2*y1 - x1*y2)
    mapping(1,2) = (u1*x2 - u2*x1) / (y1*x2 - y2*x1)

    mapping(2,1) = (v2*y1 - v1*y2) / (x2*y1 - x1*y2)
    mapping(2,2) = (v1*x2 - v2*x1) / (y1*x2 - y2*x1)

  end subroutine get_triangle_mapping!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  use_triangle_mapping
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine use_triangle_mapping(&
       u, v, &
       x, y, &
       mapping)!{{{

    real(kind=RKIND), intent(out) :: &
         u, v !< Output:

    real(kind=RKIND), intent(in) :: &
         x, y !< Input:

    real(kind=RKIND), dimension(2,2), intent(in) :: &
         mapping !< Input:

    u = mapping(1,1) * x + mapping(1,2) * y
    v = mapping(2,1) * x + mapping(2,2) * y

  end subroutine use_triangle_mapping!}}}

!-----------------------------------------------------------------------

end module cice_velocity_solver_pwl
