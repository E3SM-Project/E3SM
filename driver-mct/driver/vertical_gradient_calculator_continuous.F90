module vertical_gradient_calculator_continuous

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module defines a subclass of vertical_gradient_calculator_base_type for
  ! computing piecewise continuous vertical gradients using a matrix solve.

#include "shr_assert.h"
  use seq_comm_mct, only : logunit
  use vertical_gradient_calculator_base, only : vertical_gradient_calculator_base_type
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use mct_mod
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use shr_sys_mod, only : shr_sys_abort
  
  implicit none
  private

  public :: vertical_gradient_calculator_continuous_type

  type, extends(vertical_gradient_calculator_base_type) :: &
       vertical_gradient_calculator_continuous_type
     private

     integer :: min_elevation_class
     integer :: max_elevation_class
     integer :: nelev ! number of elevation classes
     integer :: num_points
     real(r8), allocatable :: field(:,:)   ! field(i,j) is elevation class i, field j
     real(r8), allocatable :: topo(:,:)    ! topo(i,j) is elevation class i, field j

     real(r8), allocatable :: vertical_gradient(:,:)  ! precomputed vertical gradients; vertical_gradient(i,j) is elevation class i, field j

     ! Bounds of each elevation class. This array has one more element than the number of
     ! elevation classes, since it contains lower and upper bounds for each elevation
     ! class. The indices go (min_elevation_class-1):max_elevation_class. These bounds
     ! are guaranteed to be monotonically increasing.
     real(r8), allocatable :: elevclass_bounds(:)

   contains
     procedure :: calc_vertical_gradient

     procedure, private :: set_data_from_attr_vect ! extract data from an attribute vector
     procedure, private :: precompute_vertical_gradients ! compute vertical gradients for all ECs

  end type vertical_gradient_calculator_continuous_type

  interface vertical_gradient_calculator_continuous_type
     module procedure constructor
  end interface vertical_gradient_calculator_continuous_type

contains

  !-----------------------------------------------------------------------
  function constructor(attr_vect, fieldname, toponame, &
       nelev, elevclass_names, &
       elevclass_bounds) &
       result(this)
    !
    ! !DESCRIPTION:
    ! Creates a vertical_gradient_calculator_continuous_type object by reading the
    ! necessary data from the provided attribute vector.
    !
    ! Pre-condition: elevclass_bounds must be monotonically increasing.
    !
    ! Pre-condition: Topographic heights in the attribute vector should all lie inside the
    ! bounds of their respective elevation class (given by elevclass_bounds), with the
    ! possible exception of the lowest elevation class (topographic heights can lie below
    ! the arbitrary lower bound of the elevation class) and the highest elevation class
    ! (topographic heights can lie above the arbitrary upper bound of the elevation
    ! class). For grid cells where this is not true, sets vertical gradient to 0 for all
    ! elevation classes.
    !
    ! The attribute vector is assumed to have fields named fieldname //
    ! elevclass_names(1), toponame // elevclass_names(1), etc.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(vertical_gradient_calculator_continuous_type) :: this  ! function result
    type(mct_aVect)  , intent(in) :: attr_vect           ! attribute vector in which we can find the data
    character(len=*) , intent(in) :: fieldname           ! base name of the field of interest
    character(len=*) , intent(in) :: toponame            ! base name of the topographic field
    integer          , intent(in) :: nelev               ! number of elevation classes (indexing assumed to start at 1)

    ! strings corresponding to each elevation class
    character(len=*) , intent(in) :: elevclass_names(:)

    ! bounds of each elevation class; this array should have one more element than the
    ! number of elevation classes, since it contains lower and upper bounds for each
    ! elevation class
    real(r8)         , intent(in) :: elevclass_bounds(0:)
    !
    ! !LOCAL VARIABLES:
    integer :: pt

    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(elevclass_names) == (/nelev/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(elevclass_bounds) == (/nelev/)), errMsg(__FILE__, __LINE__))

    this%nelev = nelev
    allocate(this%elevclass_bounds(0:nelev))
    this%elevclass_bounds(:) = elevclass_bounds(:)
    call this%check_elevclass_bounds_monotonic_increasing(this%elevclass_bounds)

    call this%set_data_from_attr_vect(attr_vect, fieldname, toponame, elevclass_names)

    allocate(this%vertical_gradient(this%nelev, this%num_points))
    this%vertical_gradient(:,:) = nan

    ! FIXME(wjs, 2016-04-26) Uncomment this call to check_topo - but change it so that it
    ! sets a flag, and then we'll set vertical gradients to 0 wherever topos are bad.
    !
    ! call this%check_topo()

    ! For this implementation of the vertical gradient calculator, we compute all vertical
    ! gradients in object construction. This is because we compute them all simultaneously
    ! rather than independently. (So then, the call to the routine that would normally
    ! compute vertical gradients for one elevation class simply returns the pre-computed
    ! vertical gradients for that elevation class.)

    do pt = 1, this%num_points
       call this%precompute_vertical_gradients(pt)
    end do

  end function constructor

  !-----------------------------------------------------------------------
  subroutine calc_vertical_gradient(this, elevation_class, vertical_gradient)
    !
    ! !DESCRIPTION:
    ! Returns the vertical gradient for all points, at a given elevation class.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_continuous_type), intent(in) :: this
    integer, intent(in) :: elevation_class

    ! vertical_gradient should already be allocated to the appropriate size
    real(r8), intent(out) :: vertical_gradient(:)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'calc_vertical_gradient'
    !-----------------------------------------------------------------------

    SHR_ASSERT((size(vertical_gradient) == this%num_points), errMsg(__FILE__, __LINE__))

    if (elevation_class < 1 .or. &
         elevation_class > this%nelev) then
       write(logunit,*) subname, ': ERROR: elevation class out of bounds: ', &
            elevation_class, this%nelev
       call shr_sys_abort(subname//': ERROR: elevation class out of bounds')
    end if

    vertical_gradient(:) = this%vertical_gradient(elevation_class, :)

  end subroutine calc_vertical_gradient


  !-----------------------------------------------------------------------
  subroutine set_data_from_attr_vect(this, attr_vect, fieldname, toponame, elevclass_names)
    !
    ! !DESCRIPTION:
    ! Extract data from an attribute vector.
    !
    ! Sets this%num_points, and allocates and sets this%field and this%topo.
    !
    ! TODO(wjs, 2016-04-26) The current flow is that the constructor calls this
    ! routine. It could be better to move this routine into a factory class that creates
    ! objects by (1) calling this routine to extract fields from the attribute vector, and
    ! then (2) calling the constructor of this class using these extracted data (so the
    ! constructor would never need to be passed an attribute vector).
    !
    ! !USES:
    use mct_mod
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_continuous_type), intent(inout) :: this
    type(mct_aVect)  , intent(in) :: attr_vect ! attribute vector in which we can find the data
    character(len=*) , intent(in) :: fieldname ! base name of the field of interest
    character(len=*) , intent(in) :: toponame  ! base name of the topographic field
    character(len=*) , intent(in) :: elevclass_names(:) ! strings corresponding to each elevation class
    !
    ! !LOCAL VARIABLES:
    integer :: elevclass
    character(len=:), allocatable :: fieldname_ec
    character(len=:), allocatable :: toponame_ec

    ! The following temporary array is needed because mct wants pointers
    real(r8), pointer :: temp(:)
    
    character(len=*), parameter :: subname = 'set_data_from_attr_vect'
    !-----------------------------------------------------------------------

    this%num_points = mct_aVect_lsize(attr_vect)

    allocate(this%field(this%nelev, this%num_points))
    allocate(this%topo(this%nelev, this%num_points))
    allocate(temp(this%num_points))
    
    do elevclass = 1, this%nelev
       fieldname_ec = trim(fieldname) // trim(elevclass_names(elevclass))
       call mct_aVect_exportRattr(attr_vect, fieldname_ec, temp)
       this%field(elevclass,:) = temp(:)

       toponame_ec = trim(toponame) // trim(elevclass_names(elevclass))
       call mct_aVect_exportRattr(attr_vect, toponame_ec, temp)
       this%topo(elevclass,:) = temp(:)
    end do

    deallocate(temp)
    
  end subroutine set_data_from_attr_vect

  !-----------------------------------------------------------------------
  subroutine precompute_vertical_gradients(this, pt)
    !
    ! !DESCRIPTION:
    ! Compute and save vertical gradients for all elevation classes.
    !
    ! Computes a gradient in each elevation class such that the field is continuous at
    ! interfaces and the sum over squared differences from the mean is minimized.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_continuous_type), intent(inout) :: this
    integer, intent(in) :: pt  ! point to compute gradients for (1..this%num_points)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: field(this%nelev) ! mean field value of each elevation class
    real(r8) :: topo(this%nelev)  ! mean topo of each elevation class
    real(r8) :: grad(this%nelev)  ! computed gradient
    real(r8) :: topo_interface(0:this%nelev)  ! elevations at interfaces between classes
    real(r8) :: dl(this%nelev)  ! lower 1/2 widths of elevation classes
    real(r8) :: du(this%nelev)  ! upper 1/2 widths of elevation classes
    real(r8) :: h_lo(this%nelev) ! lower bounds for computing norms
    real(r8) :: h_hi(this%nelev) ! upper bounds for computing norms
    real(r8) :: dgrad(this%nelev) ! grad - grad_mean
    real(r8) :: weight_grad(this%nelev) ! weight for dgrad in solution
    real(r8) :: grad_mean               ! mean value of gradient

    real(r8) :: diag(this%nelev-1) ! diagonal of tridiagonal matrix
    real(r8) :: subd(this%nelev-1) ! subdiagonal of tridiagonal matrix
    real(r8) :: supd(this%nelev-1) ! superdiagonal of tridiagonal matrix

    real(r8) :: b(this%nelev-1) ! rhs in A*x = b
    real(r8) :: A(this%nelev-1, this%nelev) ! matrix in A*x = b
    real(r8) :: A_AT(this%nelev-1, this%nelev-1) ! A * (transpose of A)
    real(r8) :: Tinv(this%nelev-1, this%nelev-1) ! inverse of tridiagonal matrix T = A * AT
    real(r8) :: AT_Tinv(this%nelev, this%nelev-1) ! (transpose of A) * Tinv
    real(r8) :: x_least_norm(this%nelev) ! least-norm solution x in A*x = b

    ! FIXME(wjs, 2016-04-27) Rename to nelev, or probably just delete
    integer :: n

    ! FIXME(wjs, 2016-04-26) Rename to ec
    integer :: i


    character(len=*), parameter :: subname = 'precompute_vertical_gradients'
    !-----------------------------------------------------------------------

    field(:) = this%field(:,pt)
    topo(:) = this%topo(:,pt)

    ! FIXME(wjs, 2016-04-27) Rename to elevclass_bounds? Update: just delete this
    ! temporary variable.
    topo_interface(:) = this%elevclass_bounds(:)

    n = this%nelev

    do i = 1, n
       dl(i) = topo(i) - topo_interface(i-1)  ! dl(1) is never used
       du(i) = topo_interface(i) - topo(i)    ! du(n) is never used
    end do

    ! FIXME(wjs, 2016-04-26) Extract method for this loop: returns weight_grad in each
    ! elevation class
    do i = 1, n

       if (i == 1) then
          ! If topo(1) is near the top of EC1, then the weight just includes twice the
          ! width of [topo(i) .. topo_interface(i)]
          h_lo(i) = max(topo_interface(i-1), (topo(i) - (topo_interface(i) - topo(i))))
       else
          h_lo(i) = topo_interface(i-1)
       end if

       if (i == n) then
          ! If topo(n) is near the bottom of EC N, then the weight just includes twice
          ! the width of [topo_interface(i-1) .. topo(i)]
          h_hi(i) = min(topo_interface(i), (topo(i) + (topo(i) - topo_interface(i-1))))
       else
          h_hi(i) = topo_interface(i)
       end if

       ! set gradient weights based on h_hi - h_lo in each class
       weight_grad(i) = (h_hi(i) - h_lo(i)) / (h_hi(n) - h_lo(1))

       ! FIXME(wjs, 2016-04-26) Check that weight_grad > 0.
       !       If weight_grad is just slightly > 0, the matrix will be poorly conditioned.

    end do

    !--------------------------------------------------------------------
    ! Set up matrix problem for gradient solution.
    ! The idea is to match field values at interfaces.
    !
    ! For each class n:
    !    field(n) + du(n)*grad(n) = field(n+1) - dl(n+1)*grad(n+1)
    !
    ! Rearrange to get
    !    du(n)*grad(n) + dl(n+1)*grad(n+1) = field(n+1) - field(n)
    ! 
    ! This is a bidiagonal matrix system:
    !
    !    | du(1)    dl(2)                        |   | grad(1) |    | field(2) - field(1) |
    !    |          du(2)   dl(3)                |   | grad(2) |    | field(3) - field(2) |
    !    |                  du(3)   dl(4)        | * | grad(3) | =  | field(4) - field(3) |
    !    |                          du(4)  dl(5) |   | grad(4) |    | field(5) - field(4) |
    !                                                | grad(5) |
    !
    ! The solution is underdetermined (4 equations, 5 unknowns).
    ! So we add an additional constraint: 
    ! Minimize the norm of the difference between the gradient in each class and the
    !  mean gradient, weighted by the range of the elevation class.
    ! That is, minimize the sum over i of (wt(i) * (grad(i) - grad_mean))^2.
    !
    ! The mean gradient is given by 
    !
    !                field(n) - field(1)
    !    grad_mean = __________________
    !                 topo(n) - topo(1)
    !
    ! The weights are
    !
    !         h_hi(i) - h_lo(i)
    ! wt(i) = _________________
    !         h_hi(n) - h_lo(1)
    !
    ! Putting in the weights and rearranging, we get
    !
    !         | wt(1) * dgrad(1) |    | field(2) - field(1) |       | wt(1) * grad_mean |
    !         | wt(2) * dgrad(2) |    | field(3) - field(2) |       | wt(2) * grad_mean |
    !     A * | wt(3) * dgrad(3) | =  | field(4) - field(3) | - A * | wt(3) * grad_mean |
    !         | wt(4) * dgrad(4) |    | field(5) - field(4) |       | wt(4) * grad_mean |
    !         | wt(5) * dgrad(5) |                                  | wt(5) * grad_mean |
    !
    ! where A is the bidiagonal matrix above, adjusted to include the gradient weights.
    ! E.g., du(1) becomes du(1)/wt(1), etc.
    ! This system is in the form A*x = b, where we want to solve for x.
    !
    ! It can be shown that the least-norm solution of A*x = b is given by
    !
    !    x = A^T (A*A^T)^{-1} b
    !
    ! So given A, we need to compute (A*A^T), takes its inverse, premultiply by A^T,
    ! and multiply the result by b.
    ! Given the least-norm solution, it is straightforward to compute the gradients.
    !--------------------------------------------------------------------

    ! Fill A
    ! A has (n-1) rows and n columns

    A(:,:) = 0._r8
    do i = 1, n-1
       A(i,i) = du(i) / weight_grad(i)
       A(i,i+1) = dl(i+1) / weight_grad(i+1)
    end do

    ! Compute A * A^T, a tridiagonal matrix of size (n-1)
    
    A_AT = matmul(A, transpose(A))

    ! Compute tridiagonal entries of (A * A^T)

    do i = 1, n-1
       diag(i) = A_AT(i,i)        ! diagonal
       if (i < n-1) then
          supd(i) = A_AT(i,i+1)   ! superdiagonal
          subd(i) = A_AT(i+1,i)   ! subdiagonal
       else
          supd(i) = 0._r8
          subd(i) = 0._r8
       end if
    end do

    ! Compute inverse of (A * A^T)
    ! Tinv has size (n-1,n-1)

    call tridiagonal_inverse(diag, supd, subd, Tinv)

    ! Premultiply by A^T
    ! A^T * Tinv has n rows and (n-1) columns

    AT_Tinv = matmul(transpose(A), Tinv)

    ! Compute the rhs vector b
    ! Start with the field differences

    do i = 1, n-1
       b(i) = field(i+1) - field(i)
    enddo

    ! Compute mean gradient

    grad_mean = (field(n) - field(1)) / (topo(n) - topo(1)) 

    ! Subtract (A * weight_grad*grad_mean), a vector of size(n-1)

    b(:) = b(:) - matmul(A, weight_grad(:)*grad_mean)

    ! Multiply AT_Tinv by b to get the least-norm solution
    ! b has size (n-1), x_least_norm has size n

    x_least_norm = matmul(AT_Tinv, b)

    ! Divide by the weighting factor to get dgrad
    dgrad(:) = x_least_norm(:) / weight_grad(:)

    ! Add dgrad to the mean to get the total gradient
    grad(:) = grad_mean + dgrad(:)

    this%vertical_gradient(:,pt) = grad(:)

  end subroutine precompute_vertical_gradients


end module vertical_gradient_calculator_continuous
