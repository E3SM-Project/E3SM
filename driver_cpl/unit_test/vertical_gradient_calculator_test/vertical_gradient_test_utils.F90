module vertical_gradient_test_utils

  ! Utilities to aid testing of vertical gradient calculators

#include "shr_assert.h"
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use shr_kind_mod , only : r8 => shr_kind_r8
  use avect_wrapper_mod
  use mct_mod, only : mct_aVect

  implicit none
  private

  public :: two_digit_string
  public :: create_av

contains

  function two_digit_string(val)
    ! Converts val to a two-digit string
    character(len=2) :: two_digit_string
    integer, intent(in) :: val

    write(two_digit_string, '(i2.2)') val
  end function two_digit_string

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

end module vertical_gradient_test_utils

