! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module: bitarray_module
!
! Purpose: This module provides a two-dimensional bit array and a set of 
!   routines to manipulate and examine the bits of the array.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module init_atm_bitarray

   type bitarray_type
      integer, pointer, dimension(:,:) :: iarray ! Storage array
      integer :: nx, ny                 ! Number of bits in the x and y directions
      integer :: x_int_dim, y_int_dim   ! Number of integers in the x and y directions
      integer :: integer_size           ! Number of bits in an integer
   end type bitarray_type
 
   contains
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Name: bitarray_create
   !
   ! Purpose: Allocate and initialize a bit array so that all bits are FALSE 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine bitarray_create(b, i, j)
   
      implicit none
  
      ! Arguments
      integer, intent(in) :: i, j   
      type (bitarray_type), intent(out) :: b
  
      b%integer_size = bit_size(b%integer_size) 
  
      b%nx = i
      b%ny = j
  
      b%x_int_dim = ceiling(real(b%nx)/real(b%integer_size))
      b%y_int_dim = b%ny 

      nullify(b%iarray)
      allocate(b%iarray(b%x_int_dim, b%y_int_dim))
      b%iarray = 0
 
   end subroutine bitarray_create


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Name: bitarray_copy
   !
   ! Purpose: Duplicate a bitarray.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine bitarray_copy(src, dst)
   
      implicit none
  
      ! Arguments
      type (bitarray_type), intent(in) :: src
      type (bitarray_type), intent(out) :: dst
  
      dst%integer_size = src%integer_size
  
      dst%nx = src%nx
      dst%ny = src%ny
  
      dst%x_int_dim = src%x_int_dim
      dst%y_int_dim = src%y_int_dim

      allocate(dst%iarray(dst%x_int_dim, dst%y_int_dim))
      dst%iarray = src%iarray
 
   end subroutine bitarray_copy
 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Name: bitarray_set
   !
   ! Purpose: Set the bit located at (i,j) to TRUE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine bitarray_set(b, i, j)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: i, j
      type (bitarray_type), intent(inout) :: b
  
      ! Local variables
      integer :: n_integer, n_bit
  
      n_integer = ((i-1) / b%integer_size) + 1
      n_bit = mod((i-1), b%integer_size)
  
      b%iarray(n_integer, j) = ibset(b%iarray(n_integer, j), n_bit) 
 
   end subroutine bitarray_set
 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Name: bitarray_clear
   !
   ! Purpose: Set the bit located at (i,j) to FALSE 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine bitarray_clear(b, i, j)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: i, j
      type (bitarray_type), intent(inout) :: b
  
      ! Local variables
      integer :: n_integer, n_bit
  
      n_integer = ((i-1) / b%integer_size) + 1
      n_bit = mod((i-1), b%integer_size)
  
      b%iarray(n_integer, j) = ibclr(b%iarray(n_integer, j), n_bit) 
 
   end subroutine bitarray_clear
 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Name: bitarray_test
   !
   ! Purpose: To return the value of the bit located at (i,j)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   function bitarray_test(b, i, j)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: i, j
      type (bitarray_type), intent(in) :: b
  
      ! Local variables
      logical :: bitarray_test
      integer :: n_integer, n_bit
  
      n_integer = ((i-1) / b%integer_size) + 1
      n_bit = mod((i-1), b%integer_size)
  
      bitarray_test = btest(b%iarray(n_integer,j), n_bit) 
 
   end function bitarray_test
 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Name: bitarray_merge
   !
   ! Purpose: The first bitarray argument, b1, is set to the union of the .TRUE. 
   !   bits in b1 and b2. That is, after returning, a bit x in b1 is set if
   !   either x was set in b1 or x was set in b2. Thus, b1 AND b2 MUST BE BIT 
   !   ARRAYS OF THE SAME DIMENSIONS.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine bitarray_merge(b1, b2)
 
      implicit none
  
      ! Arguments
      type (bitarray_type), intent(inout) :: b1, b2
  
      ! Local variables
      integer :: i, j
  
      if (b1%x_int_dim /= b2%x_int_dim .or. b1%y_int_dim /= b2%y_int_dim) then
         write(0,*) 'ERROR: In bitarray_merge(), b1 and b2 have different dimensions.'
      end if
  
      do i=1,b1%x_int_dim
         do j=1,b1%y_int_dim
            b1%iarray(i,j) = ior(b1%iarray(i,j), b2%iarray(i,j))  
         end do
      end do
 
   end subroutine bitarray_merge
   
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Name: bitarray_destroy
   !
   ! Purpose: To deallocate all allocated memory associated with the bit array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine bitarray_destroy(b)
 
      implicit none
  
      ! Arguments
      type (bitarray_type), intent(inout) :: b
  
      if (associated(b%iarray)) then
         deallocate(b%iarray)
      else
         write(0,*) 'WARNING: In bitarray_destroy(), b is not allocated.'
      end if
 
   end subroutine bitarray_destroy
 
end module init_atm_bitarray
