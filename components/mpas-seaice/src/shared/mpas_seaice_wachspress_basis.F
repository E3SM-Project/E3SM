!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_wachspress_basis
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 4th November 2022
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_wachspress_basis

  use mpas_derived_types

  implicit none

  private
  save

  public :: &
       seaice_calc_wachspress_coefficients, &
       seaice_wachspress_indexes, &
       seaice_wachspress_basis_function, &
       seaice_wachspress_basis_derivative

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_calc_wachspress_coefficients
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_calc_wachspress_coefficients(&
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       nCells, &
       nEdgesOnCell, &
       xLocal, &
       yLocal)!{{{

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         wachspressKappa !< Output:

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         wachspressA, & !< Output:
         wachspressB    !< Output:

    integer, intent(in) :: &
         nCells !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell !< Input:

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

  end subroutine seaice_calc_wachspress_coefficients!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_wachspress_indexes
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_wachspress_indexes(&
       nEdgesOnCell, &
       nEdgesOnCellSubset, &
       vertexIndexSubset)

    use seaice_mesh, only: &
         seaice_wrapped_index

    integer, intent(in) :: &
         nEdgesOnCell !< Input:

    integer, dimension(:), intent(out) :: &
         nEdgesOnCellSubset !< Output:

    integer, dimension(:,:), intent(out) :: &
         vertexIndexSubset !< Output:

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

  end subroutine seaice_wachspress_indexes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_wachspress_basis_function
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_wachspress_basis_function(&
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

    use seaice_mesh, only: &
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
         nEdgesOnCellSubset !< Input:

    integer, dimension(:,:), intent(in) :: &
         vertexIndexSubset !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         wachpress !< Output:

    real(kind=RKIND), dimension(size(x),nEdgesOnCell) :: &
         numerator

    real(kind=RKIND), dimension(size(x)) :: &
         denominator, &
         edgeEquation

    integer :: &
         jVertex

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

  end subroutine seaice_wachspress_basis_function!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_wachspress_basis_derivative
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_wachspress_basis_derivative(&
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
         nEdgesOnCellSubset !< Input:

    integer, dimension(:,:), intent(in) :: &
         vertexIndexSubset !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         wachspressU, & !< Output:
         wachspressV !< Output:

    real(kind=RKIND), dimension(size(x),2,nEdgesOnCell) :: &
         derivative

    real(kind=RKIND), dimension(size(x),nEdgesOnCell) :: &
         numerator

    real(kind=RKIND), dimension(size(x),2) :: &
         sum_of_derivatives, &
         sum_of_products, &
         product

    real(kind=RKIND), dimension(size(x)) :: &
         denominator, &
         edgeEquation

    integer :: &
         jVertex

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

  end subroutine seaice_wachspress_basis_derivative!}}}

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
         nEdgesOnCellSubset !< Input:

    integer, dimension(:,:), intent(in) :: &
         vertexIndexSubset !< Input:

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
         nEdgesOnCellSubset !< Input:

    integer, dimension(:,:), intent(in) :: &
         vertexIndexSubset !< Input:

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         derivative !< Output:

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         sum_of_products, & !< Input/Output:
         product            !< Input/Output:

    real(kind=RKIND), dimension(:), intent(inout) :: &
         edgeEquation !< Input/Output:

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

!-----------------------------------------------------------------------

end module seaice_wachspress_basis
