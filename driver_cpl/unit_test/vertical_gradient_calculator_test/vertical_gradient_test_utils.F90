module vertical_gradient_test_utils

  ! Utilities to aid testing of vertical gradient calculators

#include "shr_assert.h"
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use shr_kind_mod , only : r8 => shr_kind_r8
  use vertical_gradient_calculator_base, only : vertical_gradient_calculator_base_type

  implicit none
  private

  public :: all_gradients_one_point  ! Return gradients for all ECs for one point

contains

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

