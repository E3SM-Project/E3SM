module vertical_gradient_test_utils

  ! Utilities to aid testing of vertical gradient calculators

#include "shr_assert.h"
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use shr_kind_mod , only : r8 => shr_kind_r8
  use avect_wrapper_mod
  use mct_mod, only : mct_aVect
  use vertical_gradient_calculator_base, only : vertical_gradient_calculator_base_type

  implicit none
  private

  public :: two_digit_string
  public :: elevclass_names
  public :: create_av
  public :: all_gradients_one_point  ! Return gradients for all ECs for one point

contains

  function two_digit_string(val)
    ! Converts val to a two-digit string
    character(len=2) :: two_digit_string
    integer, intent(in) :: val

    write(two_digit_string, '(i2.2)') val
  end function two_digit_string

  function elevclass_names(n_elev_classes)
    ! Returns array of elevation class names
    integer, intent(in) :: n_elev_classes
    character(len=16) :: elevclass_names(n_elev_classes)

    integer :: i

    do i = 1, n_elev_classes
       elevclass_names(i) = two_digit_string(i)
    end do
  end function elevclass_names

  subroutine create_av(topo, data, toponame, dataname, av)
    ! Creates the attribute vector 'av'
    real(r8), intent(in) :: topo(:,:)  ! topo(i,j) is point i, elevation class j
    real(r8), intent(in) :: data(:,:)  ! data(i,j) is point i, elevation class j
    character(len=*), intent(in) :: toponame
    character(len=*), intent(in) :: dataname
    type(mct_aVect), intent(out) :: av

    integer :: npts
    integer :: n_elev_classes
    integer :: elevclass
    character(len=64), allocatable :: attr_tags(:)

    npts = size(topo, 1)
    n_elev_classes = size(topo, 2)

    SHR_ASSERT_ALL((ubound(data) == (/npts, n_elev_classes/)), errMsg(__FILE__, __LINE__))

    allocate(attr_tags(2*n_elev_classes))
    do elevclass = 1, n_elev_classes
       attr_tags(elevclass) = dataname // two_digit_string(elevclass)
    end do
    do elevclass = 1, n_elev_classes
       attr_tags(n_elev_classes + elevclass) = toponame // two_digit_string(elevclass)
    end do
       
    call create_aVect_with_data_rows_are_points(av, &
         attr_tags = attr_tags, &
         data = reshape([data, topo], [npts, n_elev_classes * 2]))

  end subroutine create_av

  function all_gradients_one_point(calculator, n_elev_classes, npts, pt) result(gradients)
    ! Return gradients for all ECs for one point
    class(vertical_gradient_calculator_base_type), intent(in) :: calculator
    integer, intent(in) :: n_elev_classes  ! number of elevation classes in this calculator
    integer, intent(in) :: npts  ! number of points in this calculator
    integer, intent(in) :: pt  ! point of interest
    real(r8) :: gradients(n_elev_classes)  ! function result

    integer :: ec
    real(r8) :: gradients_one_ec(npts)

    do ec = 1, n_elev_classes
       call calculator%calc_vertical_gradient(ec, gradients_one_ec)
       gradients(ec) = gradients_one_ec(pt)
    end do
  end function all_gradients_one_point

end module vertical_gradient_test_utils

