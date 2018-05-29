!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_velocity_solver_wachspress
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_velocity_solver_wachspress

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timer
  use mpas_log, only: mpas_log_write

  implicit none

  private
  save

  public :: &
       seaice_init_velocity_solver_wachspress

contains

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_velocity_solver_wachspress
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_velocity_solver_wachspress(&
       mesh, &
       velocity_variational, &
       boundary, &
       rotateCartesianGrid, &
       includeMetricTerms, &
       integrationType, &
       integrationOrder)!{{{

    use mpas_timer

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

    character(len=strKIND), intent(in) :: &
         integrationType !< Input:

    integer, intent(in) :: &
         integrationOrder !< Input:

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
         wachspressKappa, &
         integrationType, &
         integrationOrder)

    call seaice_cell_vertices_at_vertex(&
         mesh, &
         cellVerticesAtVertex)

    deallocate(xLocal)
    deallocate(yLocal)
    deallocate(wachspressKappa)
    deallocate(wachspressA)
    deallocate(wachspressB)

    call mpas_timer_stop("Velocity solver Wachpress init")

  end subroutine seaice_init_velocity_solver_wachspress!}}}

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
          wachspressA(iVertex, iCell) = &
               (yLocal(i2,iCell) - yLocal(i1,iCell)) / (xLocal(i1,iCell) * yLocal(i2,iCell) - xLocal(i2,iCell) * yLocal(i1,iCell))
          wachspressB(iVertex, iCell) = &
               (xLocal(i1,iCell) - xLocal(i2,iCell)) / (xLocal(i1,iCell) * yLocal(i2,iCell) - xLocal(i2,iCell) * yLocal(i1,iCell))

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
                  (wachspressA(i2,iCell) * (xLocal(i0,iCell) - xLocal(i1,iCell)) + &
                   wachspressB(i2,iCell) * (yLocal(i0,iCell) - yLocal(i1,iCell))) / &
                  (wachspressA(i0,iCell) * (xLocal(i1,iCell) - xLocal(i0,iCell)) + &
                   wachspressB(i0,iCell) * (yLocal(i1,iCell) - yLocal(i0,iCell)))

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

    use seaice_velocity_solver_variational_shared, only: &
         seaice_wrapped_index

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
       i2 = seaice_wrapped_index(jVertex + 1, nEdgesOnCell)

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

    use seaice_velocity_solver_variational_shared, only: &
         seaice_wrapped_index

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

    wachspressU(:) = derivative(:,1,iVertex) / denominator(:) - &
         (numerator(:,iVertex) / denominator(:)**2) * sum_of_derivatives(:,1)
    wachspressV(:) = derivative(:,2,iVertex) / denominator(:) - &
         (numerator(:,iVertex) / denominator(:)**2) * sum_of_derivatives(:,2)

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

    use seaice_velocity_solver_variational_shared, only: &
         seaice_wrapped_index

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

          basisGradientU(iBasisVertex,:,iCell) = 0.0_RKIND
          basisGradientV(iBasisVertex,:,iCell) = 0.0_RKIND

          iGradientVertex = iBasisVertex
          basisGradientU(iBasisVertex,iGradientVertex,iCell) = derivativeU(iGradientVertex)
          basisGradientV(iBasisVertex,iGradientVertex,iCell) = derivativeV(iGradientVertex)

          iGradientVertex = seaice_wrapped_index(iBasisVertex - 1, nEdgesOnCell(iCell))
          basisGradientU(iBasisVertex,iGradientVertex,iCell) = derivativeU(iGradientVertex)
          basisGradientV(iBasisVertex,iGradientVertex,iCell) = derivativeV(iGradientVertex)

          iGradientVertex = seaice_wrapped_index(iBasisVertex + 1, nEdgesOnCell(iCell))
          basisGradientU(iBasisVertex,iGradientVertex,iCell) = derivativeU(iGradientVertex)
          basisGradientV(iBasisVertex,iGradientVertex,iCell) = derivativeV(iGradientVertex)

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
       wachspressKappa, &
       integrationType, &
       integrationOrder)!{{{

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

    character(len=strKIND), intent(in) :: &
         integrationType !< Input:

    integer, intent(in) :: &
         integrationOrder !< Input:

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

    integer :: &
         nIntegrationPoints

    real(kind=RKIND), dimension(:), allocatable :: &
         integrationU, &
         integrationV, &
         integrationWeights

    real(kind=RKIND) :: &
         normalizationFactor

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)

    call get_integration_factors(&
         integrationType, &
         integrationOrder, &
         nIntegrationPoints, &
         integrationU, &
         integrationV, &
         integrationWeights, &
         normalizationFactor)

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
                  nIntegrationPoints, &
                  integrationU, &
                  integrationV, &
                  integrationWeights, &
                  normalizationFactor)

          enddo ! jVertex

       enddo ! iVertex

    enddo ! iCell

    deallocate(integrationU)
    deallocate(integrationV)
    deallocate(integrationWeights)

  end subroutine integrate_wachspress!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_integration_factors
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 18th October 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_integration_factors(integrationType, integrationOrder, nIntegrationPoints, u, v, weights, normalizationFactor)

    character(len=strKIND), intent(in) :: &
         integrationType

    integer, intent(in) :: &
         integrationOrder

    integer, intent(out) :: &
         nIntegrationPoints

    real(kind=RKIND), dimension(:), allocatable, intent(out) :: &
         u, &
         v, &
         weights

    real(kind=RKIND), intent(out) :: &
         normalizationFactor

    if (trim(integrationType) == "trapezoidal") then

       call get_integration_factors_trapezoidal(&
            integrationOrder, &
            nIntegrationPoints, &
            u, &
            v, &
            weights, &
            normalizationFactor)

    else if (trim(integrationType) == "dunavant") then

       call get_integration_factors_dunavant(&
            integrationOrder, &
            nIntegrationPoints, &
            u, &
            v, &
            weights, &
            normalizationFactor)

    else if (trim(integrationType) == "fekete") then

       call get_integration_factors_fekete(&
            integrationOrder, &
            nIntegrationPoints, &
            u, &
            v, &
            weights, &
            normalizationFactor)

    else

       ! unknown integration type
       call mpas_log_write("get_integration_factors: Unknown wachspress integration type: "//trim(integrationType), MPAS_LOG_CRIT)

    endif

  end subroutine get_integration_factors

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_integration_factors_trapezoidal
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 18th October 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_integration_factors_trapezoidal(&
       integrationOrder, &
       nIntegrationPoints, &
       u, &
       v, &
       weights, &
       normalizationFactor)

    integer, intent(in) :: &
         integrationOrder

    integer, intent(out) :: &
         nIntegrationPoints

    real(kind=RKIND), dimension(:), allocatable, intent(out) :: &
         u, &
         v, &
         weights

    real(kind=RKIND), intent(out) :: &
         normalizationFactor

    integer :: &
         nIntegrationTriangles

    integer :: &
         i, j, ij

    nIntegrationTriangles = integrationOrder

    ! total number of integration points in sub triangle
    nIntegrationPoints = ((nIntegrationTriangles+1)**2 + (nIntegrationTriangles+1)) / 2

    ! allocate integration factors
    allocate(u(nIntegrationPoints))
    allocate(v(nIntegrationPoints))
    allocate(weights(nIntegrationPoints))

    ! get integration point canonical location
    ij = 1
    do i = 0, nIntegrationTriangles
       do j = 0, nIntegrationTriangles-i

          u(ij) = real(i,RKIND) / real(nIntegrationTriangles,RKIND)
          v(ij) = real(j,RKIND) / real(nIntegrationTriangles,RKIND)

          ij = ij + 1

       enddo ! j
    enddo ! i

    ! get the weights
    ij = 1
    do i = 0, nIntegrationTriangles
       do j = 0, nIntegrationTriangles-i

          weights(ij) = 0.0_RKIND

          if (i<=nIntegrationTriangles-j) then

             if (i==nIntegrationTriangles .or. j==nIntegrationTriangles .or. (i==0 .and. j==0)) then

                weights(ij) = 1.0_RKIND

             else if ((j==0 .and. i/=0 .and. i/=nIntegrationTriangles) .or. &
                      (i==0 .and. j/=0 .and. j/=nIntegrationTriangles) .or. &
                      (i==nIntegrationTriangles-j .and. i/=0 .and. j/=0)) then

                weights(ij) = 3.0_RKIND

             else

                weights(ij) = 6.0_RKIND

             endif

          endif

          ij = ij + 1

       enddo ! j
    enddo ! i

    ! normalization factor
    normalizationFactor = 6.0_RKIND * real(nIntegrationTriangles,RKIND)**2

  end subroutine get_integration_factors_trapezoidal

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_integration_factors
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 18th October 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_integration_factors_dunavant(&
       integrationOrder, &
       nIntegrationPoints, &
       u, &
       v, &
       weights, &
       normalizationFactor)

    integer, intent(in) :: &
         integrationOrder

    integer, intent(out) :: &
         nIntegrationPoints

    real(kind=RKIND), dimension(:), allocatable, intent(out) :: &
         u, &
         v, &
         weights

    real(kind=RKIND), intent(out) :: &
         normalizationFactor

    ! D. A. Dunavant, High degree efficient symmetrical Gaussian quadrature rules for the triangle,
    ! Int. J. Num. Meth. Engng, 21(1985):1129-1148.

    normalizationFactor = 2.0_RKIND

    if (modulo(integrationOrder,2) /= 0) then
       call mpas_log_write("get_integration_factors_dunavant: odd orders of integration not recommended", MPAS_LOG_WARN)
    endif

    select case (integrationOrder)
    case(1)

       nIntegrationPoints = 1

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            0.33333333333333_RKIND /)

       v = (/ &
            0.33333333333333_RKIND /)

       weights = (/ &
            1.00000000000000_RKIND /)

    case (2)

       nIntegrationPoints = 3

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            0.16666666666667_RKIND, 0.16666666666667_RKIND, 0.66666666666667_RKIND /)

       v = (/ &
            0.16666666666667_RKIND, 0.66666666666667_RKIND, 0.16666666666667_RKIND /)

       weights = (/ &
            0.33333333333333_RKIND, 0.33333333333333_RKIND, 0.33333333333333_RKIND /)

    case (3)

       nIntegrationPoints = 4

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            0.33333333333333_RKIND, 0.20000000000000_RKIND, 0.20000000000000_RKIND, 0.60000000000000_RKIND /)

       v = (/ &
            0.33333333333333_RKIND, 0.20000000000000_RKIND, 0.60000000000000_RKIND, 0.20000000000000_RKIND /)

       weights = (/ &
           -0.56250000000000_RKIND, 0.52083333333333_RKIND, 0.52083333333333_RKIND, 0.52083333333333_RKIND /)

    case (4)

       nIntegrationPoints = 6

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            0.44594849091597_RKIND, 0.44594849091597_RKIND, 0.10810301816807_RKIND, 0.09157621350977_RKIND, &
            0.09157621350977_RKIND, 0.81684757298046_RKIND /)

       v = (/ &
            0.44594849091597_RKIND, 0.10810301816807_RKIND, 0.44594849091597_RKIND, 0.09157621350977_RKIND, &
            0.81684757298046_RKIND, 0.09157621350977_RKIND /)

       weights = (/ &
            0.22338158967801_RKIND, 0.22338158967801_RKIND, 0.22338158967801_RKIND, 0.10995174365532_RKIND, &
            0.10995174365532_RKIND, 0.10995174365532_RKIND /)

    case (5)

       nIntegrationPoints = 7

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            0.33333333333333_RKIND, 0.47014206410511_RKIND, 0.47014206410511_RKIND, 0.05971587178977_RKIND, &
            0.10128650732346_RKIND, 0.10128650732346_RKIND, 0.79742698535309_RKIND /)

       v = (/ &
            0.33333333333333_RKIND, 0.47014206410511_RKIND, 0.05971587178977_RKIND, 0.47014206410511_RKIND, &
            0.10128650732346_RKIND, 0.79742698535309_RKIND, 0.10128650732346_RKIND /)

       weights = (/ &
            0.22500000000000_RKIND, 0.13239415278851_RKIND, 0.13239415278851_RKIND, 0.13239415278851_RKIND, &
            0.12593918054483_RKIND, 0.12593918054483_RKIND, 0.12593918054483_RKIND /)

    case (6)

       nIntegrationPoints = 12

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            0.24928674517091_RKIND, 0.24928674517091_RKIND, 0.50142650965818_RKIND, 0.06308901449150_RKIND, &
            0.06308901449150_RKIND, 0.87382197101700_RKIND, 0.31035245103378_RKIND, 0.63650249912140_RKIND, &
            0.05314504984482_RKIND, 0.63650249912140_RKIND, 0.31035245103378_RKIND, 0.05314504984482_RKIND /)

       v = (/ &
            0.24928674517091_RKIND, 0.50142650965818_RKIND, 0.24928674517091_RKIND, 0.06308901449150_RKIND, &
            0.87382197101700_RKIND, 0.06308901449150_RKIND, 0.63650249912140_RKIND, 0.05314504984482_RKIND, &
            0.31035245103378_RKIND, 0.31035245103378_RKIND, 0.05314504984482_RKIND, 0.63650249912140_RKIND /)

       weights = (/ &
            0.11678627572638_RKIND, 0.11678627572638_RKIND, 0.11678627572638_RKIND, 0.05084490637021_RKIND, &
            0.05084490637021_RKIND, 0.05084490637021_RKIND, 0.08285107561837_RKIND, 0.08285107561837_RKIND, &
            0.08285107561837_RKIND, 0.08285107561837_RKIND, 0.08285107561837_RKIND, 0.08285107561837_RKIND /)

    case (7)

       nIntegrationPoints = 13

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            0.33333333333333_RKIND, 0.26034596607904_RKIND, 0.26034596607904_RKIND, 0.47930806784192_RKIND, &
            0.06513010290222_RKIND, 0.06513010290222_RKIND, 0.86973979419557_RKIND, 0.31286549600487_RKIND, &
            0.63844418856981_RKIND, 0.04869031542532_RKIND, 0.63844418856981_RKIND, 0.31286549600487_RKIND, &
            0.04869031542532_RKIND /)

       v = (/ &
            0.33333333333333_RKIND, 0.26034596607904_RKIND, 0.47930806784192_RKIND, 0.26034596607904_RKIND, &
            0.06513010290222_RKIND, 0.86973979419557_RKIND, 0.06513010290222_RKIND, 0.63844418856981_RKIND, &
            0.04869031542532_RKIND, 0.31286549600487_RKIND, 0.31286549600487_RKIND, 0.04869031542532_RKIND, &
            0.63844418856981_RKIND /)

       weights = (/ &
           -0.14957004446768_RKIND, 0.17561525743321_RKIND, 0.17561525743321_RKIND, 0.17561525743321_RKIND, &
            0.05334723560884_RKIND, 0.05334723560884_RKIND, 0.05334723560884_RKIND, 0.07711376089026_RKIND, &
            0.07711376089026_RKIND, 0.07711376089026_RKIND, 0.07711376089026_RKIND, 0.07711376089026_RKIND, &
            0.07711376089026_RKIND /)

    case (8)

       nIntegrationPoints = 16

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            0.33333333333333_RKIND, 0.45929258829272_RKIND, 0.45929258829272_RKIND, 0.08141482341455_RKIND, &
            0.17056930775176_RKIND, 0.17056930775176_RKIND, 0.65886138449648_RKIND, 0.05054722831703_RKIND, &
            0.05054722831703_RKIND, 0.89890554336594_RKIND, 0.26311282963464_RKIND, 0.72849239295540_RKIND, &
            0.00839477740996_RKIND, 0.72849239295540_RKIND, 0.26311282963464_RKIND, 0.00839477740996_RKIND /)

       v = (/ &
            0.33333333333333_RKIND, 0.45929258829272_RKIND, 0.08141482341455_RKIND, 0.45929258829272_RKIND, &
            0.17056930775176_RKIND, 0.65886138449648_RKIND, 0.17056930775176_RKIND, 0.05054722831703_RKIND, &
            0.89890554336594_RKIND, 0.05054722831703_RKIND, 0.72849239295540_RKIND, 0.00839477740996_RKIND, &
            0.26311282963464_RKIND, 0.26311282963464_RKIND, 0.00839477740996_RKIND, 0.72849239295540_RKIND /)

       weights = (/ &
            0.14431560767779_RKIND, 0.09509163426728_RKIND, 0.09509163426728_RKIND, 0.09509163426728_RKIND, &
            0.10321737053472_RKIND, 0.10321737053472_RKIND, 0.10321737053472_RKIND, 0.03245849762320_RKIND, &
            0.03245849762320_RKIND, 0.03245849762320_RKIND, 0.02723031417443_RKIND, 0.02723031417443_RKIND, &
            0.02723031417443_RKIND, 0.02723031417443_RKIND, 0.02723031417443_RKIND, 0.02723031417443_RKIND /)

    case (9)

       nIntegrationPoints = 19

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            0.333333333333333_RKIND, 0.020634961602525_RKIND, 0.489682519198738_RKIND, 0.489682519198738_RKIND, &
            0.125820817014127_RKIND, 0.437089591492937_RKIND, 0.437089591492937_RKIND, 0.623592928761935_RKIND, &
            0.188203535619033_RKIND, 0.188203535619033_RKIND, 0.910540973211095_RKIND, 0.044729513394453_RKIND, &
            0.044729513394453_RKIND, 0.036838412054736_RKIND, 0.221962989160766_RKIND, 0.036838412054736_RKIND, &
            0.741198598784498_RKIND, 0.221962989160766_RKIND, 0.741198598784498_RKIND /)

       v = (/ &
            0.333333333333333_RKIND, 0.489682519198738_RKIND, 0.020634961602525_RKIND, 0.489682519198738_RKIND, &
            0.437089591492937_RKIND, 0.125820817014127_RKIND, 0.437089591492937_RKIND, 0.188203535619033_RKIND, &
            0.623592928761935_RKIND, 0.188203535619033_RKIND, 0.044729513394453_RKIND, 0.910540973211095_RKIND, &
            0.044729513394453_RKIND, 0.221962989160766_RKIND, 0.036838412054736_RKIND, 0.741198598784498_RKIND, &
            0.036838412054736_RKIND, 0.741198598784498_RKIND, 0.221962989160766_RKIND /)

       weights = (/ &
            0.097135796282799_RKIND, 0.031334700227139_RKIND, 0.031334700227139_RKIND, 0.031334700227139_RKIND, &
            0.077827541004774_RKIND, 0.077827541004774_RKIND, 0.077827541004774_RKIND, 0.079647738927210_RKIND, &
            0.079647738927210_RKIND, 0.079647738927210_RKIND, 0.025577675658698_RKIND, 0.025577675658698_RKIND, &
            0.025577675658698_RKIND, 0.043283539377289_RKIND, 0.043283539377289_RKIND, 0.043283539377289_RKIND, &
            0.043283539377289_RKIND, 0.043283539377289_RKIND, 0.043283539377289_RKIND /)

    case (10)

       nIntegrationPoints = 25

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            0.333333333333333_RKIND, 0.028844733232685_RKIND, 0.485577633383657_RKIND, 0.485577633383657_RKIND, &
            0.781036849029926_RKIND, 0.109481575485037_RKIND, 0.109481575485037_RKIND, 0.141707219414880_RKIND, &
            0.307939838764121_RKIND, 0.141707219414880_RKIND, 0.550352941820999_RKIND, 0.307939838764121_RKIND, &
            0.550352941820999_RKIND, 0.025003534762686_RKIND, 0.246672560639903_RKIND, 0.025003534762686_RKIND, &
            0.728323904597411_RKIND, 0.246672560639903_RKIND, 0.728323904597411_RKIND, 0.009540815400299_RKIND, &
            0.066803251012200_RKIND, 0.009540815400299_RKIND, 0.923655933587500_RKIND, 0.066803251012200_RKIND, &
            0.923655933587500_RKIND /)

       v = (/ &
            0.333333333333333_RKIND, 0.485577633383657_RKIND, 0.028844733232685_RKIND, 0.485577633383657_RKIND, &
            0.109481575485037_RKIND, 0.781036849029926_RKIND, 0.109481575485037_RKIND, 0.307939838764121_RKIND, &
            0.141707219414880_RKIND, 0.550352941820999_RKIND, 0.141707219414880_RKIND, 0.550352941820999_RKIND, &
            0.307939838764121_RKIND, 0.246672560639903_RKIND, 0.025003534762686_RKIND, 0.728323904597411_RKIND, &
            0.025003534762686_RKIND, 0.728323904597411_RKIND, 0.246672560639903_RKIND, 0.066803251012200_RKIND, &
            0.009540815400299_RKIND, 0.923655933587500_RKIND, 0.009540815400299_RKIND, 0.923655933587500_RKIND, &
            0.066803251012200_RKIND /)

       weights = (/ &
            0.090817990382754_RKIND, 0.036725957756467_RKIND, 0.036725957756467_RKIND, 0.036725957756467_RKIND, &
            0.045321059435528_RKIND, 0.045321059435528_RKIND, 0.045321059435528_RKIND, 0.072757916845420_RKIND, &
            0.072757916845420_RKIND, 0.072757916845420_RKIND, 0.072757916845420_RKIND, 0.072757916845420_RKIND, &
            0.072757916845420_RKIND, 0.028327242531057_RKIND, 0.028327242531057_RKIND, 0.028327242531057_RKIND, &
            0.028327242531057_RKIND, 0.028327242531057_RKIND, 0.028327242531057_RKIND, 0.009421666963733_RKIND, &
            0.009421666963733_RKIND, 0.009421666963733_RKIND, 0.009421666963733_RKIND, 0.009421666963733_RKIND, &
            0.009421666963733_RKIND /)

    case (12)

       nIntegrationPoints = 33

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            0.023565220452390_RKIND, 0.488217389773805_RKIND, 0.488217389773805_RKIND, 0.120551215411079_RKIND, &
            0.439724392294460_RKIND, 0.439724392294460_RKIND, 0.457579229975768_RKIND, 0.271210385012116_RKIND, &
            0.271210385012116_RKIND, 0.744847708916828_RKIND, 0.127576145541586_RKIND, 0.127576145541586_RKIND, &
            0.957365299093579_RKIND, 0.021317350453210_RKIND, 0.021317350453210_RKIND, 0.115343494534698_RKIND, &
            0.275713269685514_RKIND, 0.115343494534698_RKIND, 0.608943235779788_RKIND, 0.275713269685514_RKIND, &
            0.608943235779788_RKIND, 0.022838332222257_RKIND, 0.281325580989940_RKIND, 0.022838332222257_RKIND, &
            0.695836086787803_RKIND, 0.281325580989940_RKIND, 0.695836086787803_RKIND, 0.025734050548330_RKIND, &
            0.116251915907597_RKIND, 0.025734050548330_RKIND, 0.858014033544073_RKIND, 0.116251915907597_RKIND, &
            0.858014033544073_RKIND /)

       v = (/ &
            0.488217389773805_RKIND, 0.023565220452390_RKIND, 0.488217389773805_RKIND, 0.439724392294460_RKIND, &
            0.120551215411079_RKIND, 0.439724392294460_RKIND, 0.271210385012116_RKIND, 0.457579229975768_RKIND, &
            0.271210385012116_RKIND, 0.127576145541586_RKIND, 0.744847708916828_RKIND, 0.127576145541586_RKIND, &
            0.021317350453210_RKIND, 0.957365299093579_RKIND, 0.021317350453210_RKIND, 0.275713269685514_RKIND, &
            0.115343494534698_RKIND, 0.608943235779788_RKIND, 0.115343494534698_RKIND, 0.608943235779788_RKIND, &
            0.275713269685514_RKIND, 0.281325580989940_RKIND, 0.022838332222257_RKIND, 0.695836086787803_RKIND, &
            0.022838332222257_RKIND, 0.695836086787803_RKIND, 0.281325580989940_RKIND, 0.116251915907597_RKIND, &
            0.025734050548330_RKIND, 0.858014033544073_RKIND, 0.025734050548330_RKIND, 0.858014033544073_RKIND, &
            0.116251915907597_RKIND /)

       weights = (/ &
            0.025731066440455_RKIND, 0.025731066440455_RKIND, 0.025731066440455_RKIND, 0.043692544538038_RKIND, &
            0.043692544538038_RKIND, 0.043692544538038_RKIND, 0.062858224217885_RKIND, 0.062858224217885_RKIND, &
            0.062858224217885_RKIND, 0.034796112930709_RKIND, 0.034796112930709_RKIND, 0.034796112930709_RKIND, &
            0.006166261051559_RKIND, 0.006166261051559_RKIND, 0.006166261051559_RKIND, 0.040371557766381_RKIND, &
            0.040371557766381_RKIND, 0.040371557766381_RKIND, 0.040371557766381_RKIND, 0.040371557766381_RKIND, &
            0.040371557766381_RKIND, 0.022356773202303_RKIND, 0.022356773202303_RKIND, 0.022356773202303_RKIND, &
            0.022356773202303_RKIND, 0.022356773202303_RKIND, 0.022356773202303_RKIND, 0.017316231108659_RKIND, &
            0.017316231108659_RKIND, 0.017316231108659_RKIND, 0.017316231108659_RKIND, 0.017316231108659_RKIND, &
            0.017316231108659_RKIND /)

    case default

       call mpas_log_write(&
            "get_integration_factors_dunavant: Unimplemented integration order for Dunavant wachspress integration", &
            MPAS_LOG_CRIT)

    end select

  end subroutine get_integration_factors_dunavant

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_integration_factors
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 18th October 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_integration_factors_fekete(&
       integrationOrder, &
       nIntegrationPoints, &
       u, &
       v, &
       weights, &
       normalizationFactor)

    integer, intent(in) :: &
         integrationOrder

    integer, intent(out) :: &
         nIntegrationPoints

    real(kind=RKIND), dimension(:), allocatable, intent(out) :: &
         u, &
         v, &
         weights

    real(kind=RKIND), intent(out) :: &
         normalizationFactor

    ! M. A. TAYLOR, B. A. WINGATE, AND R. E. VINCENT, (200), "AN ALGORITHM FOR COMPUTING FEKETE POINTS IN THE TRIANGLE",
    ! SIAM J. NUMER. ANAL., Vol. 38, No. 5, pp. 1707–1720

    normalizationFactor = 2.0_RKIND

    select case (integrationOrder)
    case (1)

       nIntegrationPoints = 1

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            3.33333333333333333e-01_RKIND /)

       v = (/ &
            3.33333333333333333e-01_RKIND /)

       weights = (/ &
            1.0_RKIND /)

    case (2)

       nIntegrationPoints = 3

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            1.66666666666666666e-01_RKIND, 6.66666666666666666e-01_RKIND, 1.66666666666666666e-01_RKIND /)

       v = (/ &
            6.66666666666666666e-01_RKIND, 1.66666666666666666e-01_RKIND, 1.66666666666666666e-01_RKIND /)

       weights = (/ &
            3.33333333333333333e-01_RKIND, 3.33333333333333333e-01_RKIND, 3.33333333333333333e-01_RKIND /)

    case (3:4)

       nIntegrationPoints = 6

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            9.15762135097704655e-02_RKIND, 8.16847572980458514e-01_RKIND, &
            9.15762135097710761e-02_RKIND, 1.08103018168070275e-01_RKIND, &
            4.45948490915965612e-01_RKIND, 4.45948490915964113e-01_RKIND /)

       v = (/ &
            8.16847572980458514e-01_RKIND, 9.15762135097710761e-02_RKIND, &
            9.15762135097704655e-02_RKIND, 4.45948490915964113e-01_RKIND, &
            1.08103018168070275e-01_RKIND, 4.45948490915965612e-01_RKIND /)

       weights = (/ &
            1.09951743655321843e-01_RKIND, 1.09951743655321857e-01_RKIND, &
            1.09951743655321885e-01_RKIND, 2.23381589678011389e-01_RKIND, &
            2.23381589678011527e-01_RKIND, 2.23381589678011527e-01_RKIND /)

    case (5)

       nIntegrationPoints = 10

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            0.00000000000000000e+00_RKIND, 1.00000000000000000e+00_RKIND, &
            0.00000000000000000e+00_RKIND, 2.67327353118498978e-01_RKIND, &
            6.72817552946136210e-01_RKIND, 6.49236350054349654e-02_RKIND, &
            6.71649853904175198e-01_RKIND, 6.54032456800035522e-02_RKIND, &
            2.69376706913982855e-01_RKIND, 3.38673850389605513e-01_RKIND /)

       v = (/ &
            1.00000000000000000e+00_RKIND, 0.00000000000000000e+00_RKIND, &
            0.00000000000000000e+00_RKIND, 6.72819921871012694e-01_RKIND, &
            2.67328859948191944e-01_RKIND, 6.71653011149382917e-01_RKIND, &
            6.49251690028951334e-02_RKIND, 2.69378936645285116e-01_RKIND, &
            6.54054874919145490e-02_RKIND, 3.38679989302702156e-01_RKIND /)

       weights = (/ &
            1.31356049751916795e-02_RKIND, 1.31358306034076201e-02_RKIND, &
            1.37081973800151392e-02_RKIND, 1.17419193291163376e-01_RKIND, &
            1.17420611913379477e-01_RKIND, 1.24012589655715613e-01_RKIND, &
            1.24015246126072495e-01_RKIND, 1.25930230276426303e-01_RKIND, &
            1.25933026682913923e-01_RKIND, 2.25289469095714456e-01_RKIND /)

    case (6)

       nIntegrationPoints = 11

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            5.72549866774768601e-02_RKIND, 8.95362640024579104e-01_RKIND, 6.84475748456514044e-01_RKIND, &
            6.87462559150295305e-02_RKIND, 6.15676205575839575e-01_RKIND, 6.27946141197789465e-01_RKIND, &
            6.29091383418635686e-02_RKIND, 6.83782119205099126e-02_RKIND, 2.87529458374392255e-01_RKIND, &
            3.28783556413134614e-01_RKIND, 3.12290405013644801e-01_RKIND /)

       v = (/ &
            8.95498146789879490e-01_RKIND, 6.18282212503219533e-02_RKIND, 2.33437384976827311e-02_RKIND, &
            6.00302757472630025e-02_RKIND, 3.33461808341377175e-01_RKIND, 1.59189185992151483e-01_RKIND, &
            6.55295093705452469e-01_RKIND, 3.09117685428267230e-01_RKIND, 6.36426509179620181e-01_RKIND, &
            7.70240056424634223e-02_RKIND, 3.52344786445899505e-01_RKIND /)

       weights = (/ &
            3.80680718529555623e-02_RKIND, 3.83793553077528410e-02_RKIND, 4.62004567445618367e-02_RKIND, &
            5.34675894441989999e-02_RKIND, 8.37558269657456833e-02_RKIND, 1.01644833025517037e-01_RKIND, &
            1.01861524461366940e-01_RKIND, 1.11421831660001677e-01_RKIND, 1.12009450262946106e-01_RKIND, &
            1.24787571437558295e-01_RKIND, 1.88403488837394911e-01_RKIND /)

    case (8)

       nIntegrationPoints = 16

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            7.28492392955404355e-01_RKIND, 8.39477740995753056e-03_RKIND, 2.63112829634638112e-01_RKIND, &
            8.39477740995753056e-03_RKIND, 7.28492392955404355e-01_RKIND, 2.63112829634638112e-01_RKIND, &
            5.05472283170310122e-02_RKIND, 5.05472283170309566e-02_RKIND, 8.98905543365938087e-01_RKIND, &
            4.59292588292723236e-01_RKIND, 8.14148234145536387e-02_RKIND, 4.59292588292723125e-01_RKIND, &
            1.70569307751760324e-01_RKIND, 1.70569307751760046e-01_RKIND, 6.58861384496479685e-01_RKIND, &
            3.33333333333333370e-01_RKIND /)

       v = (/ &
            8.39477740995753056e-03_RKIND, 2.63112829634638112e-01_RKIND, 7.28492392955404355e-01_RKIND, &
            7.28492392955404355e-01_RKIND, 2.63112829634638112e-01_RKIND, 8.39477740995753056e-03_RKIND, &
            5.05472283170309566e-02_RKIND, 8.98905543365938087e-01_RKIND, 5.05472283170310122e-02_RKIND, &
            8.14148234145536387e-02_RKIND, 4.59292588292723125e-01_RKIND, 4.59292588292723236e-01_RKIND, &
            6.58861384496479685e-01_RKIND, 1.70569307751760324e-01_RKIND, 1.70569307751760046e-01_RKIND, &
            3.33333333333333370e-01_RKIND /)

       weights = (/ &
            2.72303141744348991e-02_RKIND, 2.72303141744349199e-02_RKIND, 2.72303141744349199e-02_RKIND, &
            2.72303141744349789e-02_RKIND, 2.72303141744349789e-02_RKIND, 2.72303141744349997e-02_RKIND, &
            3.24584976231980793e-02_RKIND, 3.24584976231980793e-02_RKIND, 3.24584976231981001e-02_RKIND, &
            9.50916342672845638e-02_RKIND, 9.50916342672846193e-02_RKIND, 9.50916342672846193e-02_RKIND, &
            1.03217370534718286e-01_RKIND, 1.03217370534718314e-01_RKIND, 1.03217370534718314e-01_RKIND, &
            1.44315607677787283e-01_RKIND /)

    case (9)

       nIntegrationPoints = 19

       allocate(u(nIntegrationPoints))
       allocate(v(nIntegrationPoints))
       allocate(weights(nIntegrationPoints))

       u = (/ &
            2.26739052759332704e-01_RKIND, 4.77345862087794129e-02_RKIND, 2.26577168977105115e-02_RKIND, &
            9.10074385862343016e-01_RKIND, 4.41452661673673585e-02_RKIND, 4.79944340675050984e-01_RKIND, &
            7.42657808541620557e-01_RKIND, 7.43369623518591927e-01_RKIND, 2.79454959355581213e-02_RKIND, &
            3.71861932583309532e-02_RKIND, 2.22639561442096401e-01_RKIND, 1.16082059855864395e-01_RKIND, &
            4.73822270420208358e-01_RKIND, 4.77758170054016440e-01_RKIND, 6.46387881792721997e-01_RKIND, &
            2.85357695207302253e-01_RKIND, 2.04236860041029755e-01_RKIND, 1.59370884213907937e-01_RKIND, &
            3.95698265017060125e-01_RKIND /)

       v = (/ &
            0.00000000000000000e+00_RKIND, 9.16183156802148568e-01_RKIND, 7.97193825386026345e-01_RKIND, &
            4.44666861644595901e-02_RKIND, 4.81588383854628099e-02_RKIND, 5.01294615157430568e-01_RKIND, &
            3.03405081749971196e-02_RKIND, 2.22245578824042445e-01_RKIND, 5.25527023486726308e-01_RKIND, &
            2.39263537482135413e-01_RKIND, 7.29063709376736702e-01_RKIND, 6.62507673462198188e-01_RKIND, &
            4.60334709656892230e-02_RKIND, 4.01038691325781238e-01_RKIND, 1.65342747538830548e-01_RKIND, &
            4.92973630851354261e-01_RKIND, 1.19056565447230756e-01_RKIND, 3.66261159763432431e-01_RKIND, &
            2.27511600022304139e-01_RKIND /)

       weights = (/ &
            1.58676858667487208e-02_RKIND, 2.19524732703951786e-02_RKIND, 2.40354401213296598e-02_RKIND, &
            2.58522468388647786e-02_RKIND, 2.71951393759608216e-02_RKIND, 3.02097786027936584e-02_RKIND, &
            3.70093240446550606e-02_RKIND, 4.11482921825866571e-02_RKIND, 4.26331605467379776e-02_RKIND, &
            4.71413336863812371e-02_RKIND, 5.45129844125978591e-02_RKIND, 6.26632599630084636e-02_RKIND, &
            6.31379657675310846e-02_RKIND, 7.14623133641135444e-02_RKIND, 7.51048615652924606e-02_RKIND, &
            7.98259878444318866e-02_RKIND, 8.16607475819435963e-02_RKIND, 9.37481686311500140e-02_RKIND, &
            1.04838836333477403e-01_RKIND /)

    case default

       call mpas_log_write(&
            "get_integration_factors_fekete: Unimplemented integration order for Fekete wachspress integration", &
            MPAS_LOG_CRIT)

    end select

  end subroutine get_integration_factors_fekete

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
       nIntegrationPoints, &
       integrationU, &
       integrationV, &
       integrationWeights, &
       normalizationFactor)!{{{

    use seaice_constants, only: pii

    use seaice_velocity_solver_variational_shared, only: &
         seaice_wrapped_index

    real(kind=RKIND), intent(inout) :: &
         basisIntegralsU, &
         basisIntegralsV, &
         basisIntegralsMetric

    integer, intent(in) :: &
         nEdgesOnCell, &  !< Input:
         iStressVertex, & !< Input:
         iVelocityVertex  !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xLocal, & !< Input:
         yLocal    !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB    !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input:

    integer, intent(in) :: & !< Input:
         nIntegrationPoints

    real(kind=RKIND), dimension(:), intent(in) :: &
         integrationU, & !< Input:
         integrationV, & !< Input:
         integrationWeights !< Input:

    real(kind=RKIND), intent(in) :: &
         normalizationFactor !< Input:

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
         iIntegrationPoint

    real(kind=RKIND), dimension(2,2) :: &
         mapping

    real(kind=RKIND), dimension(:), allocatable :: &
         jacobian

    integer :: &
         iSubTriangle, &
         i1, &
         i2

    character(len=strKIND) :: filename

    allocate(nEdgesOnCellSubset(nEdgesOnCell))
    allocate(vertexIndexSubset(nEdgesOnCell,nEdgesOnCell))

    call wachspress_indexes(&
         nEdgesOnCell, &
         nEdgesOnCellSubset, &
         vertexIndexSubset)

    allocate(jacobian(nEdgesOnCell))

    do iSubTriangle = 1, nEdgesOnCell

       i1 = iSubTriangle
       i2 = seaice_wrapped_index(iSubTriangle + 1, nEdgesOnCell)

       call get_triangle_mapping(&
            mapping, &
            jacobian(iSubTriangle), &
            1.0_RKIND, 0.0_RKIND, &
            0.0_RKIND, 1.0_RKIND, &
            xLocal(i1), yLocal(i1), &
            xLocal(i2), yLocal(i2))

       allocate(x(nIntegrationPoints))
       allocate(y(nIntegrationPoints))

       do iIntegrationPoint = 1, nIntegrationPoints

          call use_triangle_mapping(&
               x(iIntegrationPoint),            y(iIntegrationPoint), &
               integrationU(iIntegrationPoint), integrationV(iIntegrationPoint), &
               mapping)

       enddo ! iIntegrationPoint

       allocate(stressBasisFunction(nIntegrationPoints))
       allocate(velocityBasisFunction(nIntegrationPoints))
       allocate(velocityBasisDerivativeU(nIntegrationPoints))
       allocate(velocityBasisDerivativeV(nIntegrationPoints))

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

       basisIntegralsUSubTriangle      = 0.0_RKIND
       basisIntegralsVSubTriangle      = 0.0_RKIND
       basisIntegralsMetricSubTriangle = 0.0_RKIND

       do iIntegrationPoint = 1, nIntegrationPoints

          basisIntegralsUSubTriangle = basisIntegralsUSubTriangle + &
               jacobian(iSubTriangle) * &
               integrationWeights(iIntegrationPoint) * &
               stressBasisFunction(iIntegrationPoint) * &
               velocityBasisDerivativeU(iIntegrationPoint)

          basisIntegralsVSubTriangle = basisIntegralsVSubTriangle + &
               jacobian(iSubTriangle) * &
               integrationWeights(iIntegrationPoint) * &
               stressBasisFunction(iIntegrationPoint) * &
               velocityBasisDerivativeV(iIntegrationPoint)

          basisIntegralsMetricSubTriangle = basisIntegralsMetricSubTriangle + &
               jacobian(iSubTriangle) * &
               integrationWeights(iIntegrationPoint) * &
               stressBasisFunction(iIntegrationPoint) * &
               velocityBasisFunction(iIntegrationPoint)

       enddo ! iIntegrationPoint

       basisIntegralsU      = basisIntegralsU      + basisIntegralsUSubTriangle      / normalizationFactor
       basisIntegralsV      = basisIntegralsV      + basisIntegralsVSubTriangle      / normalizationFactor
       basisIntegralsMetric = basisIntegralsMetric + basisIntegralsMetricSubTriangle / normalizationFactor

       deallocate(stressBasisFunction)
       deallocate(velocityBasisFunction)
       deallocate(velocityBasisDerivativeU)
       deallocate(velocityBasisDerivativeV)

       deallocate(x)
       deallocate(y)

    enddo ! iSubTriangle

    deallocate(nEdgesOnCellSubset)
    deallocate(vertexIndexSubset)

    deallocate(jacobian)

  end subroutine integrate_wachspress_polygon!}}}

!-----------------------------------------------------------------------

end module seaice_velocity_solver_wachspress
