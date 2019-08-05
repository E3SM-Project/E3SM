#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module sort_mod

  implicit none
  private

  public :: &
       ! sortints(a): a(r,c) is sorted along a(1,:), and a(2:r,:) are carried
       ! along. An example application is to set a(2,:) = 1:c. After sortints
       ! returns, a(2,:) contains the permutation induced by the sort of a(1,:).
       sortints

contains
  
  subroutine sortints(a)
    integer, intent(inout) :: a(:,:)

    call sortints_impl(a, size(a,1), size(a,2))
  end subroutine sortints

  subroutine sortints_impl(a, nrow, ncol)
    use iso_c_binding, only: c_int, c_size_t, c_ptr, c_sizeof, c_loc

    interface
       subroutine qsort(base, nmemb, size, cmp) bind(c)
         use iso_c_binding, only: c_ptr, c_size_t, c_int
         type (c_ptr), value, intent(in) :: base
         integer (kind=c_size_t), value, intent(in) :: nmemb, size
         interface
            function cmp(a, b) result(o) bind(c)
              use iso_c_binding, only: c_int
              integer (kind=c_int), intent(in) :: a, b
              integer (kind=c_int) :: o
            end function cmp
         end interface
       end subroutine qsort
    end interface

    integer (kind=c_int), target, intent(inout) :: a(nrow,ncol)
    integer, intent(in) :: nrow, ncol
    type (c_ptr) :: a_ptr
    integer (kind=c_size_t) :: r_s, c_s

    a_ptr = c_loc(a)
    r_s = nrow
    c_s = ncol
    call qsort(a_ptr, c_s, r_s*c_sizeof(c_int), cmp_ints)
  end subroutine sortints_impl

  function cmp_ints(a, b) result(o) bind(c)
    use iso_c_binding, only: c_int

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

end module sort_mod
