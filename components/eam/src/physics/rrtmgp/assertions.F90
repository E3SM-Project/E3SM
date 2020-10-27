module assertions
   use shr_kind_mod,   only: r8=>shr_kind_r8
   use cam_abortutils, only: endrun
   implicit none

   ! Make all routines private unless explicitly exposed as public
   private

   ! This module provides the following public routines
   public :: assert_valid, assert_range, assert

   ! Interface blocks to allow overloading procedures
   interface assert_valid
      module procedure assert_valid_0d, assert_valid_1d, assert_valid_2d, assert_valid_3d
   end interface

   interface assert_range
      module procedure assert_range_real_1d, &
                       assert_range_real_2d, &
                       assert_range_real_3d, &
                       assert_range_integer_1d, &
                       assert_range_integer_2d
   end interface

contains

   !-------------------------------------------------------------------------------
   ! assert_range checks to see if a variable values are between a specified range
   subroutine assert_range_real_1d(x, x1, x2, varname)
      real(r8), intent(in) :: x(:)
      real(r8), intent(in) :: x1, x2
      character(len=*), intent(in) :: varname
      integer :: i

      do i = 1,size(x)
         if (x(i) < x1 .or. x(i) > x2) then
            print *, varname // ' outside range at i = ', i, ': x(i) = ', x(i)
            call endrun(varname // ' outside valid range')
         end if
      end do
   end subroutine assert_range_real_1d
   !-------------------------------------------------------------------------------
   subroutine assert_range_real_2d(x, x1, x2, varname)
      real(r8), intent(in) :: x(:,:)
      real(r8), intent(in) :: x1, x2
      character(len=*), intent(in) :: varname
      integer :: i, j

      do i = 1,size(x, 1)
         do j = 1, size(x, 2)
            if (x(i,j) < x1 .or. x(i,j) > x2) then
               print *, varname // ' outside range at i,j = ', i, j, '; x(i,j) = ', x(i,j)
               call endrun(varname // ' outside valid range')
            end if
         end do
      end do
   end subroutine assert_range_real_2d
   !-------------------------------------------------------------------------------
   subroutine assert_range_real_3d(x, x1, x2, varname)
      real(r8), intent(in) :: x(:,:,:)
      real(r8), intent(in) :: x1, x2
      character(len=*), intent(in) :: varname
      integer :: i, j, k

      do i = 1,size(x, 1)
         do j = 1,size(x, 2)
            do k = 1,size(x, 3)
               if (x(i,j,k) < x1 .or. x(i,j,k) > x2) then
                  print *, varname // ' outside range at i,j,k = ', i, j, k, &
                        '; x(i,j,k) = ', x(i,j,k)
                  call endrun(varname // ' outside valid range')
               end if
            end do
         end do
      end do
   end subroutine assert_range_real_3d
   !-------------------------------------------------------------------------------
   subroutine assert_range_integer_1d(x, x1, x2, varname)
      integer, intent(in) :: x(:)
      integer, intent(in) :: x1, x2
      character(len=*), intent(in) :: varname
      integer :: i

      do i = 1,size(x)
         if (x(i) < x1 .or. x(i) > x2) then
            print *, varname // ' outside range at i = ', i, ': x(i) = ', x(i)
            call endrun(varname // ' outside valid range')
         end if
      end do
   end subroutine assert_range_integer_1d
   !-------------------------------------------------------------------------------
   subroutine assert_range_integer_2d(x, x1, x2, varname)
      integer, intent(in) :: x(:,:)
      integer, intent(in) :: x1, x2
      character(len=*), intent(in) :: varname
      integer :: i, j

      do i = 1,size(x, 1)
         do j = 1, size(x, 2)
            if (x(i,j) < x1 .or. x(i,j) > x2) then
               print *, varname // ' outside range at i,j = ', i, j, '; x(i,j) = ', x(i,j)
               call endrun(varname // ' outside valid range')
            end if
         end do
      end do
   end subroutine assert_range_integer_2d
   !-------------------------------------------------------------------------------


   !-------------------------------------------------------------------------------
   ! assert_valid routines are designed to provide a means to check that a variable
   ! does not contain either infinities or NaNs, without knowing the actual valid
   ! range.
   subroutine assert_valid_0d(x, message)
      use infnan, only: isnan, isinf

      real(r8), intent(in) :: x
      character(len=*), intent(in) :: message

      if (isnan(x)) then
         call endrun('assert_valid failed (NaN): ' // message)
      else if (isinf(x)) then
         call endrun('assert_valid failed (inf): ' // message)
      end if
   end subroutine assert_valid_0d
   !-------------------------------------------------------------------------------
   subroutine assert_valid_1d(x, message)
      use infnan, only: isnan, isinf

      real(r8), intent(in) :: x(:)
      character(len=*), intent(in) :: message
      integer :: i

      do i = 1,size(x,1)
         if (isnan(x(i))) then
            call endrun('assert_valid failed (NaN): ' // message)
         else if (isinf(x(i))) then
            call endrun('assert_valid failed (inf): ' // message)
         end if
      end do
   end subroutine assert_valid_1d
   !-------------------------------------------------------------------------------
   subroutine assert_valid_2d(x, message)
      use infnan, only: isnan, isinf

      real(r8), intent(in) :: x(:,:)
      character(len=*), intent(in) :: message
      integer :: i, j

      do i = 1,size(x,1)
         do j = 1,size(x,2)
            if (isnan(x(i,j))) then
               call endrun('assert_valid failed (NaN): ' // message)
            else if (isinf(x(i,j))) then
               call endrun('assert_valid failed (inf): ' // message)
            end if
         end do
      end do
   end subroutine assert_valid_2d
   !-------------------------------------------------------------------------------
   subroutine assert_valid_3d(x, message)
      use infnan, only: isnan, isinf

      real(r8), intent(in) :: x(:,:,:)
      character(len=*), intent(in) :: message
      integer :: i, j, k

      do i = 1,size(x,1)
         do j = 1,size(x,2)
            do k = 1,size(x,3)
               if (isnan(x(i,j,k))) then
                  call endrun('assert_valid failed (NaN): ' // message)
               else if (isinf(x(i,j,k))) then
                  call endrun('assert_valid failed (inf): ' // message)
               end if
            end do
         end do
      end do
   end subroutine assert_valid_3d
   !-------------------------------------------------------------------------------


   ! Assert just checks that condition is true, and aborts execution of the
   ! program if it is not.
   subroutine assert(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         call endrun('Assertion failed: ' // message)
      end if
   end subroutine
   !-------------------------------------------------------------------------------

end module assertions
