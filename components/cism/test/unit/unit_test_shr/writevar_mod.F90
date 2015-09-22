! This module contains subroutines for writing variables in a compiler and machine-
! independent way, so that we can diff output that has been created with different
! compilers / machines.

module writevar_mod
   use shr_kind_mod

   implicit none
   save

   private


   ! public interfaces
   public :: writevar
   
   interface writevar
      module procedure writevar_char
      module procedure writevar_int
      module procedure writevar_real
      module procedure writevar_double
      module procedure writevar_logical
   end interface writevar

contains

   subroutine writevar_char(var, name, prefix, unit)
      character(len=*), intent(in) :: var    ! variable whose value we are outputting
      character(len=*), intent(in) :: name   ! name of the variable
      character(len=*), intent(in) :: prefix ! prefix to write at the start of the line
                                             ! (we do NOT trim this)
      integer         , intent(in) :: unit   ! output unit

      write(unit, '(a, a, " = ", a)') prefix, trim(name), trim(var)
   end subroutine writevar_char

   subroutine writevar_int(var, name, prefix, unit)
      integer         , intent(in) :: var    ! variable whose value we are outputting
      character(len=*), intent(in) :: name   ! name of the variable
      character(len=*), intent(in) :: prefix ! prefix to write at the start of the line
                                             ! (we do NOT trim this)
      integer         , intent(in) :: unit   ! output unit

      write(unit, '(a, a, " = ", i0)') prefix, trim(name), var
   end subroutine writevar_int

   subroutine writevar_real(var, name, prefix, unit)
      real(SHR_KIND_R4), intent(in) :: var   ! variable whose value we are outputting
      character(len=*), intent(in) :: name   ! name of the variable
      character(len=*), intent(in) :: prefix ! prefix to write at the start of the line
                                             ! (we do NOT trim this)
      integer         , intent(in) :: unit   ! output unit

      write(unit, '(a, a, " = ", f0.6)') prefix, trim(name), var
   end subroutine writevar_real

   subroutine writevar_double(var, name, prefix, unit)
      real(SHR_KIND_R8), intent(in) :: var   ! variable whose value we are outputting
      character(len=*), intent(in) :: name   ! name of the variable
      character(len=*), intent(in) :: prefix ! prefix to write at the start of the line
                                             ! (we do NOT trim this)
      integer         , intent(in) :: unit   ! output unit

      write(unit, '(a, a, " = ", f0.12)') prefix, trim(name), var
   end subroutine writevar_double

   subroutine writevar_logical(var, name, prefix, unit)
      logical         , intent(in) :: var    ! variable whose value we are outputting
      character(len=*), intent(in) :: name   ! name of the variable
      character(len=*), intent(in) :: prefix ! prefix to write at the start of the line
                                             ! (we do NOT trim this)
      integer         , intent(in) :: unit   ! output unit

      write(unit, '(a, a, " = ", l7)') prefix, trim(name), var
   end subroutine writevar_logical
end module writevar_mod
