! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2019,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
module mo_util_array
!
! This module provide utilites for sanitizing input arrays:
!    checking values and sizes
! These are in a module so code can be written for both CPUs and GPUs
! Used only by Fortran classes so routines don't need C bindings and can use assumed-shape
! Currently only for 3D arrays; could extend through overloading to other ranks
!
!
  use mo_rte_kind,      only: wp
  implicit none
  private
  public :: any_vals_less_than, any_vals_outside
contains
!-------------------------------------------------------------------------------------------------
  logical function any_vals_less_than(array, minVal)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp),                   intent(in) :: minVal

    real(wp) :: minValue
    integer  :: i, j, k

    ! This could be written far more compactly as
    !       any_vals_less_than = any(array < minVal)
    ! but an explicit loop also works on GPUs
    minValue = minVal
    !$acc parallel loop collapse(3) copyin(array) reduction(min:minValue)
    do k = 1, size(array,3)
      do j = 1, size(array,2)
        do i = 1, size(array,1)
          minValue = min(array(i,j,k), minValue)
        end do
      end do
    end do
    any_vals_less_than = (minValue < minVal)
  end function any_vals_less_than
  ! ---------------------------------
  logical function any_vals_outside(array, minVal, maxVal)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp),                   intent(in) :: minVal, maxVal

    real(wp) :: minValue, maxValue
    integer  :: i, j, k

    ! This could be written far more compactly as
    !   any_vals_outside = any(array < minVal .or. array > maxVal)
    ! but an explicit loop also works on GPUs
    minValue = minVal
    maxValue = maxVal
    !$acc parallel loop collapse(3) copyin(array) reduction(min:minValue) reduction(max:maxValue)
    do k = 1, size(array,3)
      do j = 1, size(array,2)
        do i = 1, size(array,1)
          minValue = min(array(i,j,k), minValue)
          maxValue = max(array(i,j,k), minValue)
        end do
      end do
    end do
    any_vals_outside = (minValue < minVal .or. maxValue > maxVal)
  end function any_vals_outside
! ---------------------------------
end module mo_util_array
