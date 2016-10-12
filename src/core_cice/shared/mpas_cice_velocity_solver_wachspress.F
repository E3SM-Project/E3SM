!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_velocity_solver_wachspress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_velocity_solver_wachspress

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timer

  implicit none

  private
  save

  public :: &
       cice_init_velocity_solver_wachspress

contains

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_velocity_solver_wachspress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_init_velocity_solver_wachspress(&
       mesh, &
       velocity_variational, &
       boundary, &
       rotateCartesianGrid, &
       includeMetricTerms)!{{{

    use mpas_timer

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

    integer :: &
         iCell, &
         iVertex

    integer, pointer :: &
         nCells, &
         maxEdges, i1, i2

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex

    real(kind=RKIND), dimension(:), pointer :: &
         tanLatVertexRotatedOverRadius

    real(kind=RKIND), dimension(:,:), allocatable :: &
         xLocal, &
         yLocal

    real(kind=RKIND), dimension(:,:), allocatable :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND), dimension(:,:,:), allocatable :: &
         wachspressKappa

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisGradientU, &
         basisGradientV, &
         basisIntegralsU, &
         basisIntegralsV, &
         basisIntegralsMetric

    call mpas_timer_start("Velocity solver Wachpress init")

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "maxEdges", maxEdges)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)

    call MPAS_pool_get_array(velocity_variational, "cellVerticesAtVertex", cellVerticesAtVertex)
    call MPAS_pool_get_array(velocity_variational, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)
    call MPAS_pool_get_array(velocity_variational, "basisGradientU", basisGradientU)
    call MPAS_pool_get_array(velocity_variational, "basisGradientV", basisGradientV)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsU", basisIntegralsU)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsV", basisIntegralsV)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsMetric", basisIntegralsMetric)

    allocate(xLocal(maxEdges,nCells))
    allocate(yLocal(maxEdges,nCells))
    allocate(wachspressKappa(maxEdges,maxEdges,nCells))
    allocate(wachspressA(maxEdges,nCells))
    allocate(wachspressB(maxEdges,nCells))

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

    call calc_wachspress_coefficients(&
         mesh, &
         wachspressKappa, &
         wachspressA, &
         wachspressB, &
         xLocal, &
         yLocal)

    call calculate_wachspress_derivatives(&
         mesh, &
         basisGradientU, &
         basisGradientV, &
         xLocal, &
         yLocal, &
         wachspressA, &
         wachspressB, &
         wachspressKappa)

    call integrate_wachspress(&
         mesh, &
         basisIntegralsU, &
         basisIntegralsV, &
         basisIntegralsMetric, &
         xLocal, &
         yLocal, &
         wachspressA, &
         wachspressB, &
         wachspressKappa)

    call cice_cell_vertices_at_vertex(&
         mesh, &
         cellVerticesAtVertex)

    deallocate(xLocal)
    deallocate(yLocal)
    deallocate(wachspressKappa)
    deallocate(wachspressA)
    deallocate(wachspressB)

    call mpas_timer_stop("Velocity solver Wachpress init")

  end subroutine cice_init_velocity_solver_wachspress!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calc_wachspress_coefficients
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine calc_wachspress_coefficients(&
       mesh, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       xLocal, &
       yLocal)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         wachspressKappa !< Output:

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         wachspressA, & !< Output:
         wachspressB    !< Output:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         xLocal, & !< Input:
         yLocal    !< Input:

    integer :: &
         iCell, &
         iVertex, &
         i0, &
         i1, &
         i2, &
         jVertex

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)

    ! loop over cells
    do iCell = 1, nCells

       ! loop over vertices
       do iVertex = 1, nEdgesOnCell(iCell)

          ! end points of line segment
          i1 = iVertex - 1
          i2 = iVertex
          if (i1 < 1) i1 = i1 + nEdgesOnCell(iCell)

          ! solve for the line segment equation
          wachspressA(iVertex, iCell) = (yLocal(i2,iCell) - yLocal(i1,iCell)) / (xLocal(i1,iCell) * yLocal(i2,iCell) - xLocal(i2,iCell) * yLocal(i1,iCell))
          wachspressB(iVertex, iCell) = (xLocal(i1,iCell) - xLocal(i2,iCell)) / (xLocal(i1,iCell) * yLocal(i2,iCell) - xLocal(i2,iCell) * yLocal(i1,iCell))

       enddo ! iVertex

       ! loop over vertices
       do iVertex = 1, nEdgesOnCell(iCell)

          ! determine kappa
          wachspressKappa(1,iVertex,iCell) = 1.0_RKIND

          do jVertex = 2, nEdgesOnCell(iCell)

             ! previous, this and next vertex
             i0 = jVertex - 1
             i1 = jVertex
             i2 = jVertex + 1
             if (i2 > nEdgesOnCell(iCell)) i2 = i2 - nEdgesOnCell(iCell)

             wachspressKappa(jVertex,iVertex,iCell) = wachspressKappa(jVertex-1,iVertex,iCell) * &
                  (wachspressA(i2,iCell) * (xLocal(i0,iCell) - xLocal(i1,iCell)) + wachspressB(i2,iCell) * (yLocal(i0,iCell) - yLocal(i1,iCell))) / &
                  (wachspressA(i0,iCell) * (xLocal(i1,iCell) - xLocal(i0,iCell)) + wachspressB(i0,iCell) * (yLocal(i1,iCell) - yLocal(i0,iCell)))

          enddo ! jVertex

       enddo ! iVertex

    enddo ! iCell

  end subroutine calc_wachspress_coefficients!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wachspress_indexes
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine wachspress_indexes(&
       nEdgesOnCell, &
       nEdgesOnCellSubset, &
       vertexIndexSubset)

    use cice_velocity_solver_variational_shared, only: &
         cice_wrapped_index

    integer, intent(in) :: &
         nEdgesOnCell

    integer, dimension(:), intent(out) :: &
         nEdgesOnCellSubset

    integer, dimension(:,:), intent(out) :: &
         vertexIndexSubset

    integer :: &
         jVertex, &
         kVertex, &
         i1, i2

    do jVertex = 1, nEdgesOnCell

       i1 = jVertex
       i2 = cice_wrapped_index(jVertex + 1, nEdgesOnCell)

       nEdgesOnCellSubset(jVertex) = 0

       do kVertex = 1, nEdgesOnCell

          if (kVertex /= i1 .and. kVertex /= i2) then
             nEdgesOnCellSubset(jVertex) = nEdgesOnCellSubset(jVertex) + 1
             vertexIndexSubset(jVertex,nEdgesOnCellSubset(jVertex)) = kVertex
          endif

       enddo ! kVertex

    enddo ! jVertex

  end subroutine wachspress_indexes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wachspress_basis_function
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine wachspress_basis_function(&
       nEdgesOnCell, &
       iVertex, &
       x, &
       y, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       nEdgesOnCellSubset, &
       vertexIndexSubset, &
       wachpress)!{{{

    use cice_velocity_solver_variational_shared, only: &
         cice_wrapped_index

    integer, intent(in) :: &
         nEdgesOnCell, & !< Input:
         iVertex         !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         x, & !< Input:
         y    !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB    !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCellSubset

    integer, dimension(:,:), intent(in) :: &
         vertexIndexSubset

    real(kind=RKIND), dimension(:), intent(out) :: &
         wachpress !< Output:

    real(kind=RKIND), dimension(:,:), allocatable :: &
         numerator

    real(kind=RKIND), dimension(:), allocatable :: &
         denominator, &
         edgeEquation

    integer :: &
         jVertex

    allocate(denominator(size(x)))
    allocate(numerator(size(x),nEdgesOnCell))
    allocate(edgeEquation(size(x)))

    ! sum over numerators to get denominator
    denominator(:) = 0.0_RKIND

    do jVertex = 1, nEdgesOnCell

      call wachspress_numerator(&
           nEdgesOnCell, jVertex, iVertex, x(:), y(:), &
           wachspressKappa, wachspressA, wachspressB, &
           nEdgesOnCellSubset, vertexIndexSubset, &
           edgeEquation(:), &
           numerator(:,jVertex))

       denominator(:) = denominator(:) + numerator(:,jVertex)

    enddo ! jVertex

    wachpress(:) = numerator(:,iVertex) / denominator(:)

    deallocate(denominator)
    deallocate(numerator)
    deallocate(edgeEquation)

  end subroutine wachspress_basis_function!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wachspress_basis_derivative
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine wachspress_basis_derivative(&
       nEdgesOnCell, &
       iVertex, &
       x, &
       y, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       nEdgesOnCellSubset, &
       vertexIndexSubset, &
       wachspressU, &
       wachspressV)!{{{

    integer, intent(in) :: &
         nEdgesOnCell, & !< Input:
         iVertex         !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         x, & !< Input:
         y    !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB    !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCellSubset

    integer, dimension(:,:), intent(in) :: &
         vertexIndexSubset

    real(kind=RKIND), dimension(:), intent(out) :: &
         wachspressU, & !< Output:
         wachspressV !< Output:

    real(kind=RKIND), dimension(:,:,:), allocatable :: &
         derivative

    real(kind=RKIND), dimension(:,:), allocatable :: &
         numerator, &
         sum_of_derivatives, &
         sum_of_products, &
         product

    real(kind=RKIND), dimension(:), allocatable :: &
         denominator, &
         edgeEquation

    integer :: &
         jVertex

    allocate(denominator(size(x)))
    allocate(sum_of_derivatives(size(x),2))
    allocate(numerator(size(x),nEdgesOnCell))
    allocate(derivative(size(x),2,nEdgesOnCell))
    allocate(sum_of_products(size(x),2))
    allocate(product(size(x),2))
    allocate(edgeEquation(size(x)))

    ! sum over numerators to get denominator
    denominator(:) = 0.0_RKIND
    sum_of_derivatives(:,:) = 0.0_RKIND

    do jVertex = 1, nEdgesOnCell

       call wachspress_numerator(&
            nEdgesOnCell, jVertex, iVertex, x(:), y(:), &
            wachspressKappa, wachspressA, wachspressB, &
            nEdgesOnCellSubset, vertexIndexSubset, &
            edgeEquation, &
            numerator(:,jVertex))

       denominator(:) = denominator(:) + numerator(:,jVertex)

       call wachspress_numerator_derivative(&
            nEdgesOnCell, jVertex, iVertex, x(:), y(:), &
            wachspressKappa, wachspressA, wachspressB, &
            nEdgesOnCellSubset, vertexIndexSubset, &
            sum_of_products, product, edgeEquation, &
            derivative(:,:,jVertex))

       sum_of_derivatives(:,:) = sum_of_derivatives(:,:) + derivative(:,:,jVertex)

    enddo ! jVertex

    wachspressU(:) = derivative(:,1,iVertex) / denominator(:) - (numerator(:,iVertex) / denominator(:)**2) * sum_of_derivatives(:,1)
    wachspressV(:) = derivative(:,2,iVertex) / denominator(:) - (numerator(:,iVertex) / denominator(:)**2) * sum_of_derivatives(:,2)

    deallocate(denominator)
    deallocate(sum_of_derivatives)
    deallocate(numerator)
    deallocate(derivative)
    deallocate(sum_of_products)
    deallocate(product)
    deallocate(edgeEquation)

  end subroutine wachspress_basis_derivative!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wachspress_numerator
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine wachspress_numerator(&
       nEdgesOnCell, &
       jVertex, &
       iVertex, &
       x, &
       y, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       nEdgesOnCellSubset, &
       vertexIndexSubset, &
       edgeEquation, &
       numerator)!{{{

    integer, intent(in) :: &
         nEdgesOnCell, & !< Input:
         jVertex, &      !< Input:
         iVertex         !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         x, & !< Input:
         y    !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB    !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCellSubset

    integer, dimension(:,:), intent(in) :: &
         vertexIndexSubset

    real(kind=RKIND), dimension(:), intent(inout) :: &
         edgeEquation

    real(kind=RKIND), dimension(:), intent(out) :: &
         numerator !< Output:

    integer :: &
         kVertex

    numerator(:) = 1.0_RKIND

    do kVertex = 1, nEdgesOnCellSubset(jVertex)

       call wachspress_edge_equation(&
            x(:), y(:), &
            wachspressA(vertexIndexSubset(jVertex,kVertex)), &
            wachspressB(vertexIndexSubset(jVertex,kVertex)), &
            edgeEquation(:))

       numerator(:) = numerator(:) * edgeEquation(:)

    enddo ! jVertex

    numerator(:) = numerator(:) * wachspressKappa(jVertex,iVertex)

  end subroutine wachspress_numerator!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wachspress_numerator_derivative
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine wachspress_numerator_derivative(&
       nEdgesOnCell, &
       jVertex, &
       iVertex, &
       x, &
       y, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       nEdgesOnCellSubset, &
       vertexIndexSubset, &
       sum_of_products, &
       product, &
       edgeEquation, &
       derivative)!{{{

    integer, intent(in) :: &
         nEdgesOnCell, & !< Input:
         jVertex, &      !< Input:
         iVertex         !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         x, & !< Input:
         y    !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB    !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCellSubset

    integer, dimension(:,:), intent(in) :: &
         vertexIndexSubset

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         derivative !< Output:

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         sum_of_products, &
         product

    real(kind=RKIND), dimension(:), intent(inout) :: &
         edgeEquation

    integer :: &
         kVertex, &
         lVertex

    sum_of_products(:,:) = 0.0_RKIND

    do kVertex = 1, nEdgesOnCellSubset(jVertex)

       product(:,:) = 1.0_RKIND

       ! lVertex < kVertex
       do lVertex = 1, kVertex - 1

          call wachspress_edge_equation(&
               x(:), y(:), &
               wachspressA(vertexIndexSubset(jVertex,lVertex)), &
               wachspressB(vertexIndexSubset(jVertex,lVertex)), &
               edgeEquation(:))

          product(:,1) = product(:,1) * edgeEquation(:)
          product(:,2) = product(:,2) * edgeEquation(:)

       enddo ! lVertex

       ! lVertex == kVertex
       product(:,1) = product(:,1) * (-wachspressA(vertexIndexSubset(jVertex,kVertex)))
       product(:,2) = product(:,2) * (-wachspressB(vertexIndexSubset(jVertex,kVertex)))

       ! lVertex > kVertex
       do lVertex = kVertex + 1, nEdgesOnCellSubset(jVertex)

          call wachspress_edge_equation(&
               x(:), y(:), &
               wachspressA(vertexIndexSubset(jVertex,lVertex)), &
               wachspressB(vertexIndexSubset(jVertex,lVertex)), &
               edgeEquation(:))

          product(:,1) = product(:,1) * edgeEquation(:)
          product(:,2) = product(:,2) * edgeEquation(:)

       enddo ! lVertex

       sum_of_products(:,:) = sum_of_products(:,:) + product(:,:)

    enddo ! jVertex

    derivative(:,:) = sum_of_products(:,:) * wachspressKappa(jVertex,iVertex)

  end subroutine wachspress_numerator_derivative!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wachspress_edge_equation
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine wachspress_edge_equation(&
       x, &
       y, &
       wachspressA, &
       wachspressB, &
       edgeEquation)

    real(kind=RKIND), dimension(:), intent(in) :: &
         x, & !< Input:
         y    !< Input:

    real(kind=RKIND), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB    !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         edgeEquation !< Output:

    edgeEquation(:) = 1.0_RKIND - wachspressA * x(:) - wachspressB * y(:)

  end subroutine wachspress_edge_equation!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calculate_wachspress_derivatives
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine calculate_wachspress_derivatives(&
       mesh, &
       basisGradientU, &
       basisGradientV, &
       xLocal, &
       yLocal, &
       wachspressA, &
       wachspressB, &
       wachspressKappa)!{{{

    use mpas_timer

    ! basisGradientUV(jVertexOnCell,iVertexOnCell,iCell)
    ! iCell         : The cell the gradients are based in
    ! iVertexOnCell : The vertex basis function the gradient is calculated from
    ! jVertexOnCell : The vertex location the gradients are calculated at

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         basisGradientU, & !< Output:
         basisGradientV    !< Output:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB, & !< Input:
         xLocal,      & !< Input:
         yLocal         !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         wachspressKappa !< Input:

    integer :: &
         iCell, &
         iBasisVertex, &
         iGradientVertex

    integer, pointer :: &
         nCells, &
         maxEdges

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    integer, dimension(:), allocatable :: &
         nEdgesOnCellSubset

    integer, dimension(:,:), allocatable :: &
         vertexIndexSubset

    real(kind=RKIND), dimension(:), allocatable :: &
         x, y, derivativeU, derivativeV

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "maxEdges", maxEdges)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)

    allocate(x(maxEdges))
    allocate(y(maxEdges))

    allocate(derivativeU(maxEdges))
    allocate(derivativeV(maxEdges))

    allocate(nEdgesOnCellSubset(maxEdges))
    allocate(vertexIndexSubset(maxEdges,maxEdges))

    ! loop over cells
    do iCell = 1, nCells

       call wachspress_indexes(&
            nEdgesOnCell(iCell), &
            nEdgesOnCellSubset(1:nEdgesOnCell(iCell)), &
            vertexIndexSubset(1:nEdgesOnCell(iCell),1:nEdgesOnCell(iCell)))

       ! loop over vertices again - derivative position
       do iGradientVertex = 1, nEdgesOnCell(iCell)

          x(iGradientVertex) = xLocal(iGradientVertex,iCell)
          y(iGradientVertex) = yLocal(iGradientVertex,iCell)

       enddo ! iGradientVertex

       ! loop over vertices - basis function
       do iBasisVertex = 1, nEdgesOnCell(iCell)

          call wachspress_basis_derivative(&
               nEdgesOnCell(iCell), &
               iBasisVertex, &
               x(1:nEdgesOnCell(iCell)), &
               y(1:nEdgesOnCell(iCell)), &
               wachspressKappa(:,:,iCell), &
               wachspressA(:,iCell), &
               wachspressB(:,iCell), &
               nEdgesOnCellSubset(1:nEdgesOnCell(iCell)), &
               vertexIndexSubset(1:nEdgesOnCell(iCell),1:nEdgesOnCell(iCell)), &
               derivativeU(1:nEdgesOnCell(iCell)), &
               derivativeV(1:nEdgesOnCell(iCell)))

          ! loop over vertices again - derivative position
          do iGradientVertex = 1, nEdgesOnCell(iCell)

             basisGradientU(iBasisVertex,iGradientVertex,iCell) = derivativeU(iGradientVertex)
             basisGradientV(iBasisVertex,iGradientVertex,iCell) = derivativeV(iGradientVertex)

          enddo ! iGradientVertex

       enddo ! iBasisVertex

    enddo ! iCell

    deallocate(nEdgesOnCellSubset)
    deallocate(vertexIndexSubset)

    deallocate(x)
    deallocate(y)

    deallocate(derivativeU)
    deallocate(derivativeV)

  end subroutine calculate_wachspress_derivatives!}}}

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
       jacobian, &
       x1, y1, &
       x2, y2, &
       u1, v1, &
       u2, v2)!{{{

    real(kind=RKIND), dimension(2,2), intent(out) :: &
         mapping !< Output:

    real(kind=RKIND), intent(out) :: &
         jacobian !< Output:

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

    jacobian = mapping(1,1) * mapping(2,2) - mapping(1,2) * mapping(2,1)

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

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  integrate_wachspress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine integrate_wachspress(&
       mesh, &
       basisIntegralsU, &
       basisIntegralsV, &
       basisIntegralsMetric, &
       xLocal, &
       yLocal, &
       wachspressA, &
       wachspressB, &
       wachspressKappa)!{{{

    use mpas_timer

    ! basisIntegralsUV (iStressVertex,iVelocityVertex,iCell)
    ! iCell         : cell integrals are performed on
    ! iStressVertex : vertex number of Wachspress function
    ! iVelocityVertex : vertex number of Wachspress derivative function
    ! Sij

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         basisIntegralsU, &   !< Output:
         basisIntegralsV, &   !< Output:
         basisIntegralsMetric !< Output:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         xLocal, &      !< Input:
         yLocal, &      !< Input:
         wachspressA, & !< Input:
         wachspressB    !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         wachspressKappa !< Input:

    real(kind=RKIND) :: &
         integration

    integer :: &
         iCell, &
         iStressVertex, &
         iVelocityVertex

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, parameter :: &
         nIntegrationTriangles = 5

    real(kind=RKIND), dimension(0:nIntegrationTriangles,0:nIntegrationTriangles) :: &
         integrationPointScaling

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)

    call setIntegrationPointScaling(&
         integrationPointScaling, &
         nIntegrationTriangles)

    do iCell = 1, nCells

       do iStressVertex = 1, nEdgesOnCell(iCell)

          do iVelocityVertex = 1, nEdgesOnCell(iCell)

             basisIntegralsU(iStressVertex,iVelocityVertex,iCell)      = 0.0_RKIND
             basisIntegralsV(iStressVertex,iVelocityVertex,iCell)      = 0.0_RKIND
             basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) = 0.0_RKIND

             call integrate_wachspress_polygon(&
                  basisIntegralsU(iStressVertex,iVelocityVertex,iCell), &
                  basisIntegralsV(iStressVertex,iVelocityVertex,iCell), &
                  basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell), &
                  nEdgesOnCell(iCell), &
                  iStressVertex, &
                  iVelocityVertex, &
                  xLocal(:,iCell), &
                  yLocal(:,iCell), &
                  wachspressA(:,iCell), &
                  wachspressB(:,iCell), &
                  wachspressKappa(:,:,iCell), &
                  nIntegrationTriangles, &
                  integrationPointScaling)

          enddo ! jVertex

       enddo ! iVertex

    enddo ! iCell

  end subroutine integrate_wachspress!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  setIntegrationPointScaling
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date December 16th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine setIntegrationPointScaling(&
       integrationPointScaling, &
       nIntegrationTriangles)

    real(kind=RKIND), dimension(0:,0:), intent(inout) :: &
         integrationPointScaling

    integer, intent(in) :: &
         nIntegrationTriangles

    integer :: &
         i, j

    do i = 0, nIntegrationTriangles
       do j = 0, nIntegrationTriangles

          integrationPointScaling(i,j) = 0.0_RKIND

          if (i<=nIntegrationTriangles-j) then

             if (i==nIntegrationTriangles .or. j==nIntegrationTriangles .or. (i==0 .and. j==0)) then

                integrationPointScaling(i,j) = 1.0_RKIND

             else if ((j==0 .and. i/=0 .and. i/=nIntegrationTriangles) .or. &
                      (i==0 .and. j/=0 .and. j/=nIntegrationTriangles) .or. &
                      (i==nIntegrationTriangles-j .and. i/=0 .and. j/=0)) then

                integrationPointScaling(i,j) = 3.0_RKIND

             else

                integrationPointScaling(i,j) = 6.0_RKIND

             endif

          endif

       enddo ! j
    enddo ! i

  end subroutine setIntegrationPointScaling

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  integrate_wachspress_polygon
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine integrate_wachspress_polygon(&
       basisIntegralsU, &
       basisIntegralsV, &
       basisIntegralsMetric, &
       nEdgesOnCell, &
       iStressVertex, &
       iVelocityVertex, &
       xLocal, &
       yLocal, &
       wachspressA, &
       wachspressB, &
       wachspressKappa, &
       nIntegrationTriangles, &
       integrationPointScaling)!{{{

    use cice_velocity_solver_variational_shared, only: &
         cice_wrapped_index

    real(kind=RKIND), intent(inout) :: &
         basisIntegralsU, &
         basisIntegralsV, &
         basisIntegralsMetric

    integer, intent(in) :: &
         nEdgesOnCell, &  !< Input:
         iStressVertex, & !< Input:
         iVelocityVertex, & !< Input:
         nIntegrationTriangles

    real(kind=RKIND), dimension(:), intent(in) :: &
         xLocal, & !< Input:
         yLocal    !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB    !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input:

    real(kind=RKIND), dimension(0:,0:), intent(in) :: &
         integrationPointScaling

    integer, dimension(:), allocatable :: &
         nEdgesOnCellSubset

    integer, dimension(:,:), allocatable :: &
         vertexIndexSubset

    real(kind=RKIND) :: &
         basisIntegralsUSubTriangle, &
         basisIntegralsVSubTriangle, &
         basisIntegralsMetricSubTriangle

    real(kind=RKIND), dimension(:), allocatable :: &
         x, &
         y, &
         stressBasisFunction, &
         velocityBasisFunction, &
         velocityBasisDerivativeU, &
         velocityBasisDerivativeV

    real(kind=RKIND) :: &
         u, &
         v

    integer :: &
         i, j, ij, nij

    real(kind=RKIND), dimension(2,2) :: &
         mapping

    real(kind=RKIND), dimension(:), allocatable :: &
         jacobian

    integer :: &
         iSubTriangle, &
         i1, &
         i2

    allocate(nEdgesOnCellSubset(nEdgesOnCell))
    allocate(vertexIndexSubset(nEdgesOnCell,nEdgesOnCell))

    call wachspress_indexes(&
         nEdgesOnCell, &
         nEdgesOnCellSubset, &
         vertexIndexSubset)

    allocate(jacobian(nEdgesOnCell))

    nij = ((nIntegrationTriangles+1)**2 + (nIntegrationTriangles+1)) / 2
    nij = nij * nEdgesOnCell
    allocate(x(nij))
    allocate(y(nij))

    ij = 1
    do iSubTriangle = 1, nEdgesOnCell

       i1 = iSubTriangle
       i2 = cice_wrapped_index(iSubTriangle + 1, nEdgesOnCell)

       call get_triangle_mapping(&
            mapping, &
            jacobian(iSubTriangle), &
            1.0_RKIND, 0.0_RKIND, &
            0.0_RKIND, 1.0_RKIND, &
            xLocal(i1), yLocal(i1), &
            xLocal(i2), yLocal(i2))

       do i = 0, nIntegrationTriangles
          do j = 0, nIntegrationTriangles-i

             u = real(i,RKIND) / real(nIntegrationTriangles,RKIND)
             v = real(j,RKIND) / real(nIntegrationTriangles,RKIND)

             call use_triangle_mapping(&
                  x(ij), y(ij), &
                  u, v, &
                  mapping)

             ij = ij + 1

          enddo ! i
       enddo ! j

    enddo ! iSubTriangle

    allocate(stressBasisFunction(nij))
    allocate(velocityBasisFunction(nij))
    allocate(velocityBasisDerivativeU(nij))
    allocate(velocityBasisDerivativeV(nij))

    call wachspress_basis_function(&
         nEdgesOnCell, iStressVertex, x, y, &
         wachspressKappa, wachspressA, wachspressB, &
         nEdgesOnCellSubset, vertexIndexSubset, &
         stressBasisFunction)

    call wachspress_basis_function(&
         nEdgesOnCell, iVelocityVertex, x, y, &
         wachspressKappa, wachspressA, wachspressB, &
         nEdgesOnCellSubset, vertexIndexSubset, &
         velocityBasisFunction)

    call wachspress_basis_derivative(&
         nEdgesOnCell, iVelocityVertex, x, y, &
         wachspressKappa, wachspressA, wachspressB, &
         nEdgesOnCellSubset, vertexIndexSubset, &
         velocityBasisDerivativeU, &
         velocityBasisDerivativeV)

    ij = 1
    do iSubTriangle = 1, nEdgesOnCell

       basisIntegralsUSubTriangle      = 0.0_RKIND
       basisIntegralsVSubTriangle      = 0.0_RKIND
       basisIntegralsMetricSubTriangle = 0.0_RKIND

       do i = 0, nIntegrationTriangles
          do j = 0, nIntegrationTriangles-i

             basisIntegralsUSubTriangle = basisIntegralsUSubTriangle + &
                  jacobian(iSubTriangle) * &
                  integrationPointScaling(i,j) * &
                  stressBasisFunction(ij) * &
                  velocityBasisDerivativeU(ij)

             basisIntegralsVSubTriangle = basisIntegralsVSubTriangle + &
                  jacobian(iSubTriangle) * &
                  integrationPointScaling(i,j) * &
                  stressBasisFunction(ij) * &
                  velocityBasisDerivativeV(ij)

             basisIntegralsMetricSubTriangle = basisIntegralsMetricSubTriangle + &
                  jacobian(iSubTriangle) * &
                  integrationPointScaling(i,j) * &
                  stressBasisFunction(ij) * &
                  velocityBasisFunction(ij)

             ij = ij + 1

          enddo ! i
       enddo ! j

       basisIntegralsU      = basisIntegralsU      + basisIntegralsUSubTriangle      / (6.0_RKIND * real(nIntegrationTriangles,RKIND)**2)
       basisIntegralsV      = basisIntegralsV      + basisIntegralsVSubTriangle      / (6.0_RKIND * real(nIntegrationTriangles,RKIND)**2)
       basisIntegralsMetric = basisIntegralsMetric + basisIntegralsMetricSubTriangle / (6.0_RKIND * real(nIntegrationTriangles,RKIND)**2)

    enddo ! iSubTriangle

    deallocate(stressBasisFunction)
    deallocate(velocityBasisFunction)
    deallocate(velocityBasisDerivativeU)
    deallocate(velocityBasisDerivativeV)

    deallocate(x)
    deallocate(y)

    deallocate(nEdgesOnCellSubset)
    deallocate(vertexIndexSubset)

    deallocate(jacobian)

  end subroutine integrate_wachspress_polygon!}}}

!-----------------------------------------------------------------------
! Hex specific plotting stuff
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_wachpress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine plot_wachpress(&
       mesh, &
       iCell, &
       nEdgesOnCell, &
       iVertex, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       xLocal, &
       yLocal)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    integer, intent(in) :: &
         iCell, &        !< Input:
         nEdgesOnCell, & !< Input:
         iVertex         !< Input:

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         wachspressKappa !< Output:

    real(kind=RKIND), dimension(:), intent(out) :: &
         wachspressA, & !< Output:
         wachspressB    !< Output:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xLocal, & !< Input:
         yLocal    !< Input:

    real(kind=RKIND) :: &
         xmin, xmax, ymin, ymax

    real(kind=RKIND) :: &
         jacobian

    real(kind=RKIND), dimension(2,2) :: &
         mapping

    integer :: &
         iSubTriangle, &
         i1, &
         i2, &
         iVertexOnCell, &
         iVertex2, &
         iObject

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

    iObject = 1

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

    open(55,file="wachspress.txt")

    write(55,fmt='(a,f10.2,a,f10.2,a)') "set xrange [",xmin,":",xmax,"]"
    write(55,fmt='(a,f10.2,a,f10.2,a)') "set yrange [",ymin,":",ymax,"]"

    ! loop over subtriangles
    do iSubTriangle = 1, nEdgesOnCell

       i1 = iSubTriangle
       i2 = iSubTriangle + 1
       if (i2 > nEdgesOnCell) i2 = i2 - nEdgesOnCell

       ! get the triangle mapping
       call get_triangle_mapping(mapping, &
                                 jacobian, &
                                 1.0_RKIND,0.0_RKIND,&
                                 0.0_RKIND,1.0_RKIND,&
                                 xLocal(i1),yLocal(i1),&
                                 xLocal(i2),yLocal(i2))

       ! plot from subtriangle
       call plot_subtriangle_2(nEdgesOnCell, iSubTriangle, iVertex, &
                             wachspressKappa, wachspressA, wachspressB, &
                             mapping, xCell(iCell), yCell(iCell), iObject)

    enddo ! iSubTriangle

    close(55)

  end subroutine plot_wachpress!}}}


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_subtriangle
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine plot_subtriangle(&
       nEdgesOnCell, &
       iSubTriangle, &
       iVertex, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       mapping, &
       x0, &
       y0, &
       iObject)!{{{

    integer, intent(in) :: &
         nEdgesOnCell, & !< Input:
         iSubTriangle, & !< Input:
         iVertex         !< Input:

    integer, intent(inout) :: &
         iObject         !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB    !< Input:

    real(kind=RKIND), dimension(2,2), intent(in) :: &
         mapping !< Input:

    real(kind=RKIND), intent(in) :: &
         x0, & !< Input:
         y0    !< Input:

    real(kind=RKIND) :: &
         u, v, &
         x, y, &
         wachspress, &
         wachspress1

    integer :: &
         i, j

    real(kind=RKIND) :: &
         x1, x2, x3, x4, x5, &
         y1, y2, y3, y4, y5, &
         d, &
         dx, dy

    integer, dimension(:), allocatable :: &
         nEdgesOnCellSubset

    integer, dimension(:,:), allocatable :: &
         vertexIndexSubset

    integer, parameter :: n = 5

    logical, parameter :: &
         lfinitedifference = .false.

    real(kind=RKIND), dimension(1) :: &
         xbasis, ybasis, var

    allocate(nEdgesOnCellSubset(nEdgesOnCell))
    allocate(vertexIndexSubset(nEdgesOnCell,nEdgesOnCell))

    call wachspress_indexes(&
         nEdgesOnCell, &
         nEdgesOnCellSubset, &
         vertexIndexSubset)

    d = 1.0_RKIND / real(n, RKIND)

    do i = 0, n
       do j = 0, n

          if (i<=n-j-1) then

             iObject = iObject + 1

             u = real(i,RKIND) / real(n,RKIND)
             v = real(j,RKIND) / real(n,RKIND)

             call use_triangle_mapping(x, y, u, v, mapping)

             if (lfinitedifference) then

                xbasis(1) = x ; ybasis(1) = y
                call wachspress_basis_function(&
                     nEdgesOnCell, iVertex, xbasis, ybasis, &
                     wachspressKappa, wachspressA, wachspressB, &
                     nEdgesOnCellSubset, vertexIndexSubset, var)
                wachspress = var(1)

                dx = 0.0_RKIND
                dy = 100.0_RKIND

                x1 = x + dx
                y1 = y + dy

                xbasis(1) = x1 ; ybasis(1) = y1
                call wachspress_basis_function(&
                     nEdgesOnCell, iVertex, xbasis, ybasis, &
                     wachspressKappa, wachspressA, wachspressB, &
                     nEdgesOnCellSubset, vertexIndexSubset, var)
                wachspress1 = var(1)

                wachspress = (wachspress1 - wachspress) / dy

             else

                xbasis(1) = x ; ybasis(1) = y
                call wachspress_basis_function(&
                     nEdgesOnCell, iVertex, xbasis, ybasis, &
                     wachspressKappa, wachspressA, wachspressB, &
                     nEdgesOnCellSubset, vertexIndexSubset, var)
                wachspress = var(1)

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

                write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') "set object ",iObject," polygon from ",&
                     x1,",",y1," to ",x2,",",y2," to ",x3,",",y3," to ",x4,",",y4," to ",x5,",",y5

                write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', wachspress, ' fillstyle solid'

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

                write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') "set object ",iObject," polygon from ",&
                     x1,",",y1," to ",x2,",",y2," to ",x4,",",y4," to ",x5,",",y5

                write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', wachspress, ' fillstyle solid'

             endif

          endif

       enddo ! j
    enddo ! i

    deallocate(nEdgesOnCellSubset)
    deallocate(vertexIndexSubset)

  end subroutine plot_subtriangle!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_subtriangle_2
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 31 October 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine plot_subtriangle_2(&
       nEdgesOnCell, &
       iSubTriangle, &
       iVertex, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       mapping, &
       x0, &
       y0, &
       iObject)!{{{

    integer, intent(in) :: &
         nEdgesOnCell, & !< Input:
         iSubTriangle, & !< Input:
         iVertex         !< Input:

    integer, intent(inout) :: &
         iObject         !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB    !< Input:

    real(kind=RKIND), dimension(2,2), intent(in) :: &
         mapping !< Input:

    real(kind=RKIND), intent(in) :: &
         x0, & !< Input:
         y0    !< Input:

    integer :: &
         i, j, &
         ii, jj

    real(kind=RKIND) :: &
         u, v, &
         uu0, vv0, &
         uu1, vv1, &
         uu2, vv2, &
         xx0, yy0, &
         xx1, yy1, &
         xx2, yy2, &
         f0, f1, f2, &
         f_avg, &
         d, &
         u0, u1, u2, &
         v0, v1, v2, &
         uuu, vvv, &
         x, y

    integer, dimension(:), allocatable :: &
         nEdgesOnCellSubset

    integer, dimension(:,:), allocatable :: &
         vertexIndexSubset

    integer, parameter :: n = 50
    integer, parameter :: nn = 1

    logical, parameter :: &
         lfinitedifference = .false.

    real(kind=RKIND), dimension(1) :: &
         xbasis, ybasis, var

    d = (1.0_RKIND / real(n,RKIND))

    allocate(nEdgesOnCellSubset(nEdgesOnCell))
    allocate(vertexIndexSubset(nEdgesOnCell,nEdgesOnCell))

    call wachspress_indexes(&
         nEdgesOnCell, &
         nEdgesOnCellSubset, &
         vertexIndexSubset)

    ! loop over integration segments
    do i = 0, n
       do j = 0, n

          ! forward integration triangle
          if (i+j<=n-1) then

             ! subtriangle starting point
             u0 = real(i,RKIND) / real(n,RKIND)
             v0 = real(j,RKIND) / real(n,RKIND)
             call use_triangle_mapping(x, y, u0, v0, mapping)

             xbasis(1) = x ; ybasis(1) = y
             call wachspress_basis_function(&
                  nEdgesOnCell, iVertex, xbasis, ybasis, &
                  wachspressKappa, wachspressA, wachspressB, &
                  nEdgesOnCellSubset, vertexIndexSubset, var)
             f0 = var(1)

             u1 = real(i+1,RKIND) / real(n,RKIND)
             v1 = real(j  ,RKIND) / real(n,RKIND)
             call use_triangle_mapping(x, y, u1, v1, mapping)

             xbasis(1) = x ; ybasis(1) = y
             call wachspress_basis_function(&
                  nEdgesOnCell, iVertex, xbasis, ybasis, &
                  wachspressKappa, wachspressA, wachspressB, &
                  nEdgesOnCellSubset, vertexIndexSubset, var)
             f1 = var(1)

             u2 = real(i  ,RKIND) / real(n,RKIND)
             v2 = real(j+1,RKIND) / real(n,RKIND)
             call use_triangle_mapping(x, y, u2, v2, mapping)

             xbasis(1) = x ; ybasis(1) = y
             call wachspress_basis_function(&
                  nEdgesOnCell, iVertex, xbasis, ybasis, &
                  wachspressKappa, wachspressA, wachspressB, &
                  nEdgesOnCellSubset, vertexIndexSubset, var)
             f2 = var(1)

             u = u0
             v = v0

             ! break up the forward triangle for sub-plotting
             do ii = 0, nn
                do jj = 0, nn

                   ! forward triangle
                   if (ii+jj<=nn-1) then

                      uu0 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv0 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu1 = real(ii+1,RKIND) / real(nn,RKIND)
                      vv1 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu2 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv2 = real(jj+1,RKIND) / real(nn,RKIND)

                      ! subtriangle location in canonical coordinates
                      uuu = (uu0 + uu1 + uu2) / 3.0_RKIND
                      vvv = (vv0 + vv1 + vv2) / 3.0_RKIND

                      ! linear basis
                      f_avg = f0 * (-uuu-vvv+1) + &
                              f1 * uuu + &
                              f2 * vvv

                      uu0 = u + uu0 * d
                      vv0 = v + vv0 * d
                      call use_triangle_mapping(xx0, yy0, uu0, vv0, mapping)

                      uu1 = u + uu1 * d
                      vv1 = v + vv1 * d
                      call use_triangle_mapping(xx1, yy1, uu1, vv1, mapping)

                      uu2 = u + uu2 * d
                      vv2 = v + vv2 * d
                      call use_triangle_mapping(xx2, yy2, uu2, vv2, mapping)

                      xx0 = x0 + xx0 ; yy0 = y0 + yy0
                      xx1 = x0 + xx1 ; yy1 = y0 + yy1
                      xx2 = x0 + xx2 ; yy2 = y0 + yy2

                      write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') "set object ",iObject," polygon from ",&
                           xx0,",",yy0," to ",xx1,",",yy1," to ",xx2,",",yy2," to ",xx0,",",yy0

                      write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', f_avg, ' fillstyle solid'

                      iObject = iObject + 1

                   endif

                   ! backward triangle
                   if (ii+jj<=nn .and. ii/=0 .and. jj/=0) then

                      uu0 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv0 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu1 = real(ii-1,RKIND) / real(nn,RKIND)
                      vv1 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu2 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv2 = real(jj-1,RKIND) / real(nn,RKIND)

                      ! subtriangle location in canonical coordinates
                      uuu = (uu0 + uu1 + uu2) / 3.0_RKIND
                      vvv = (vv0 + vv1 + vv2) / 3.0_RKIND

                      ! linear basis
                      f_avg = f0 * (-uuu-vvv+1) + &
                              f1 * uuu + &
                              f2 * vvv

                      uu0 = u + uu0 * d
                      vv0 = v + vv0 * d
                      call use_triangle_mapping(xx0, yy0, uu0, vv0, mapping)

                      uu1 = u + uu1 * d
                      vv1 = v + vv1 * d
                      call use_triangle_mapping(xx1, yy1, uu1, vv1, mapping)

                      uu2 = u + uu2 * d
                      vv2 = v + vv2 * d
                      call use_triangle_mapping(xx2, yy2, uu2, vv2, mapping)

                      xx0 = x0 + xx0 ; yy0 = y0 + yy0
                      xx1 = x0 + xx1 ; yy1 = y0 + yy1
                      xx2 = x0 + xx2 ; yy2 = y0 + yy2

                      write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') "set object ",iObject," polygon from ",&
                           xx0,",",yy0," to ",xx1,",",yy1," to ",xx2,",",yy2," to ",xx0,",",yy0

                      write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', f_avg, ' fillstyle solid'

                      iObject = iObject + 1

                   endif

                enddo ! jj
             enddo ! ii

          endif

          ! backward integration triangle
          if (i+j<=n .and. i/=0 .and. j/=0) then

             ! subtriangle starting point
             u0 = real(i,RKIND) / real(n,RKIND)
             v0 = real(j,RKIND) / real(n,RKIND)
             call use_triangle_mapping(x, y, u0, v0, mapping)

             xbasis(1) = x ; ybasis(1) = y
             call wachspress_basis_function(&
                  nEdgesOnCell, iVertex, xbasis, ybasis, &
                  wachspressKappa, wachspressA, wachspressB, &
                  nEdgesOnCellSubset, vertexIndexSubset, var)
             f0 = var(1)

             u1 = real(i-1,RKIND) / real(n,RKIND)
             v1 = real(j  ,RKIND) / real(n,RKIND)
             call use_triangle_mapping(x, y, u1, v1, mapping)

             xbasis(1) = x ; ybasis(1) = y
             call wachspress_basis_function(&
                  nEdgesOnCell, iVertex, xbasis, ybasis, &
                  wachspressKappa, wachspressA, wachspressB, &
                  nEdgesOnCellSubset, vertexIndexSubset, var)
             f1 = var(1)

             u2 = real(i  ,RKIND) / real(n,RKIND)
             v2 = real(j-1,RKIND) / real(n,RKIND)
             call use_triangle_mapping(x, y, u2, v2, mapping)

             xbasis(1) = x ; ybasis(1) = y
             call wachspress_basis_function(&
                  nEdgesOnCell, iVertex, xbasis, ybasis, &
                  wachspressKappa, wachspressA, wachspressB, &
                  nEdgesOnCellSubset, vertexIndexSubset, var)
             f2 = var(1)

             u = u0 - d
             v = v0 - d

             ! break up the backward triangle for sub-plotting
             do ii = 0, nn
                do jj = 0, nn

                   ! forward triangle
                   if (ii+jj>=nn .and. ii/=nn .and. jj/=nn) then

                      uu0 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv0 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu1 = real(ii+1,RKIND) / real(nn,RKIND)
                      vv1 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu2 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv2 = real(jj+1,RKIND) / real(nn,RKIND)

                      ! subtriangle location in canonical coordinates
                      uuu = (uu0 + uu1 + uu2) / 3.0_RKIND
                      vvv = (vv0 + vv1 + vv2) / 3.0_RKIND

                      uuu = 1.0_RKIND - uuu
                      vvv = 1.0_RKIND - vvv

                      ! linear basis
                      f_avg = f0 * (-uuu-vvv+1) + &
                              f1 * uuu + &
                              f2 * vvv

                      uu0 = u + uu0 * d
                      vv0 = v + vv0 * d
                      call use_triangle_mapping(xx0, yy0, uu0, vv0, mapping)

                      uu1 = u + uu1 * d
                      vv1 = v + vv1 * d
                      call use_triangle_mapping(xx1, yy1, uu1, vv1, mapping)

                      uu2 = u + uu2 * d
                      vv2 = v + vv2 * d
                      call use_triangle_mapping(xx2, yy2, uu2, vv2, mapping)

                      xx0 = x0 + xx0 ; yy0 = y0 + yy0
                      xx1 = x0 + xx1 ; yy1 = y0 + yy1
                      xx2 = x0 + xx2 ; yy2 = y0 + yy2

                      write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') "set object ",iObject," polygon from ",&
                           xx0,",",yy0," to ",xx1,",",yy1," to ",xx2,",",yy2," to ",xx0,",",yy0

                      write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', f_avg, ' fillstyle solid'

                      iObject = iObject + 1

                   endif

                   ! backward triangle
                   if (ii+jj>=nn+1) then

                      uu0 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv0 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu1 = real(ii-1,RKIND) / real(nn,RKIND)
                      vv1 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu2 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv2 = real(jj-1,RKIND) / real(nn,RKIND)

                      ! subtriangle location in canonical coordinates
                      uuu = (uu0 + uu1 + uu2) / 3.0_RKIND
                      vvv = (vv0 + vv1 + vv2) / 3.0_RKIND

                      uuu = 1.0_RKIND - uuu
                      vvv = 1.0_RKIND - vvv

                      ! linear basis
                      f_avg = f0 * (-uuu-vvv+1) + &
                              f1 * uuu + &
                              f2 * vvv

                      uu0 = u + uu0 * d
                      vv0 = v + vv0 * d
                      call use_triangle_mapping(xx0, yy0, uu0, vv0, mapping)

                      uu1 = u + uu1 * d
                      vv1 = v + vv1 * d
                      call use_triangle_mapping(xx1, yy1, uu1, vv1, mapping)

                      uu2 = u + uu2 * d
                      vv2 = v + vv2 * d
                      call use_triangle_mapping(xx2, yy2, uu2, vv2, mapping)

                      xx0 = x0 + xx0 ; yy0 = y0 + yy0
                      xx1 = x0 + xx1 ; yy1 = y0 + yy1
                      xx2 = x0 + xx2 ; yy2 = y0 + yy2

                      write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') "set object ",iObject," polygon from ",&
                           xx0,",",yy0," to ",xx1,",",yy1," to ",xx2,",",yy2," to ",xx0,",",yy0

                      write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', f_avg, ' fillstyle solid'

                      iObject = iObject + 1

                   endif

                enddo ! jj
             enddo ! ii

          endif

       enddo ! j
    enddo ! i

    deallocate(nEdgesOnCellSubset)
    deallocate(vertexIndexSubset)

  end subroutine plot_subtriangle_2!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_wachspress_square
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine plot_wachspress_square(&
       mesh, &
       wachspressKappa, &
       wachspressA, &
       wachspressB)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         wachspressKappa !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB    !< Input:

    integer :: &
         iCell, &
         iVertex, &
         iVertexOnCell

    real(kind=RKIND) :: &
         x, y, x0, y0, dx, dy, wachspress

    integer :: i, j

    integer, parameter :: ni = 100
    integer, parameter :: nj = 100

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex

    integer, dimension(:), allocatable :: &
         nEdgesOnCellSubset

    integer, dimension(:,:), allocatable :: &
         vertexIndexSubset

    real(kind=RKIND), dimension(1) :: &
         xbasis, ybasis

    real(kind=RKIND), dimension(1) :: &
         varU, varV

    ! init variables
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nedgesOnCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)

    iCell = 1

    x0 = 1e30
    y0 = 1e30

    do iVertexOnCell = 1, nEdgesOnCell(iCell)

       iVertex = verticesOnCell(iVertexOnCell,iCell)

       x0 = min(x0,xVertex(iVertex))
       y0 = min(y0,yVertex(iVertex))

    enddo !

    dx = 16000.0_RKIND / real(ni,RKIND)
    dy = 16000.0_RKIND / real(nj,RKIND)

    iVertexOnCell = 1

    allocate(nEdgesOnCellSubset(nEdgesOnCell(iCell)))
    allocate(vertexIndexSubset(nEdgesOnCell(iCell),nEdgesOnCell(iCell)))

    call wachspress_indexes(&
         nEdgesOnCell(iCell), &
         nEdgesOnCellSubset, &
         vertexIndexSubset)

    do i = 1, ni
       do j = 1, nj

          x = x0 + (real(i,RKIND) - 0.5_RKIND) * dx - 24000.0_RKIND
          y = x0 + (real(j,RKIND) - 0.5_RKIND) * dy - 24000.0_RKIND

          xbasis(1) = x ; ybasis(1) = y
          call wachspress_basis_derivative(&
               nEdgesOnCell(iCell), &
               iVertexOnCell, &
               xbasis, &
               ybasis, &
               wachspressKappa(:,:,iCell), &
               wachspressA(:,iCell), &
               wachspressB(:,iCell), &
               nEdgesOnCellSubset, &
               vertexIndexSubset, &
               varU, &
               varV)
          wachspress = varU(1)

          write(55,*) x+0.5_RKIND*dx,x-0.5_RKIND*dx,y+0.5_RKIND*dy,y-0.5_RKIND*dy, wachspress

       enddo ! j
    enddo ! i

    deallocate(nEdgesOnCellSubset)
    deallocate(vertexIndexSubset)

  end subroutine plot_wachspress_square!}}}

!-----------------------------------------------------------------------

end module cice_velocity_solver_wachspress
