#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module sort_mod

  use iso_c_binding, only: c_int, c_long_long, c_size_t, c_ptr, c_funptr, &
       c_sizeof, c_loc, c_funloc

  implicit none
  private

  public :: &
       ! sortints(a): a(r,c) is sorted along a(1,:), and a(2:r,:) are carried
       ! along. An example application is to set a(2,:) = 1:c. After sortints
       ! returns, a(2,:) contains the permutation induced by the sort of a(1,:).
       ! If r >= 2, ties in a(1,:) are broken by a(2,:), so the result does not
       ! depend on the qsort implementation when a(2,:) holds unique values.
       ! Generic over default (c_int) and 8-byte (c_long_long) integer arrays.
       sortints

  interface sortints
     module procedure sortints_i4, sortints_i8
  end interface sortints

  ! Single C binding for qsort; the comparator is passed as a C function
  ! pointer so the same interface serves both key widths.
  interface
     subroutine qsort(base, nmemb, size, cmp) bind(c, name='qsort')
       import :: c_ptr, c_size_t, c_funptr
       type (c_ptr), value, intent(in) :: base
       integer (kind=c_size_t), value, intent(in) :: nmemb, size
       type (c_funptr), value, intent(in) :: cmp
     end subroutine qsort
  end interface

contains

  subroutine sortints_i4(a)
    integer (kind=c_int), target, intent(inout) :: a(:,:)

    call sortints_i4_impl(a, size(a,1), size(a,2))
  end subroutine sortints_i4

  subroutine sortints_i4_impl(a, nrow, ncol)
    integer, intent(in) :: nrow, ncol
    integer (kind=c_int), target, intent(inout) :: a(nrow,ncol)

    if (nrow >= 2) then
       call qsort(c_loc(a), int(ncol, c_size_t), &
            int(nrow, c_size_t)*c_sizeof(a(1,1)), c_funloc(cmp_ints_tb))
    else
       call qsort(c_loc(a), int(ncol, c_size_t), &
            int(nrow, c_size_t)*c_sizeof(a(1,1)), c_funloc(cmp_ints))
    end if
  end subroutine sortints_i4_impl

  subroutine sortints_i8(a)
    integer (kind=c_long_long), target, intent(inout) :: a(:,:)

    call sortints_i8_impl(a, size(a,1), size(a,2))
  end subroutine sortints_i8

  subroutine sortints_i8_impl(a, nrow, ncol)
    integer, intent(in) :: nrow, ncol
    integer (kind=c_long_long), target, intent(inout) :: a(nrow,ncol)

    if (nrow >= 2) then
       call qsort(c_loc(a), int(ncol, c_size_t), &
            int(nrow, c_size_t)*c_sizeof(a(1,1)), c_funloc(cmp_ints_i8_tb))
    else
       call qsort(c_loc(a), int(ncol, c_size_t), &
            int(nrow, c_size_t)*c_sizeof(a(1,1)), c_funloc(cmp_ints_i8))
    end if
  end subroutine sortints_i8_impl

  function cmp_ints(a, b) result(o) bind(c)
    integer (kind=c_int), intent(in) :: a, b
    integer (kind=c_int) :: o

    if (a < b) then
       o = -1
    elseif (a == b) then
       o = 0
    else
       o = 1
    end if
  end function cmp_ints

  function cmp_ints_i8(a, b) result(o) bind(c)
    integer (kind=c_long_long), intent(in) :: a, b
    integer (kind=c_int) :: o

    if (a < b) then
       o = -1
    elseif (a == b) then
       o = 0
    else
       o = 1
    end if
  end function cmp_ints_i8

  ! Comparators for nrow >= 2: qsort passes the address of a column of a, so
  ! a(1:2) is (key, tie-breaker). Compare the key, then the tie-breaker.
  function cmp_ints_tb(a, b) result(o) bind(c)
    integer (kind=c_int), intent(in) :: a(2), b(2)
    integer (kind=c_int) :: o

    if (a(1) < b(1)) then
       o = -1
    elseif (a(1) > b(1)) then
       o = 1
    elseif (a(2) < b(2)) then
       o = -1
    elseif (a(2) > b(2)) then
       o = 1
    else
       o = 0
    end if
  end function cmp_ints_tb

  function cmp_ints_i8_tb(a, b) result(o) bind(c)
    integer (kind=c_long_long), intent(in) :: a(2), b(2)
    integer (kind=c_int) :: o

    if (a(1) < b(1)) then
       o = -1
    elseif (a(1) > b(1)) then
       o = 1
    elseif (a(2) < b(2)) then
       o = -1
    elseif (a(2) > b(2)) then
       o = 1
    else
       o = 0
    end if
  end function cmp_ints_i8_tb

end module sort_mod
