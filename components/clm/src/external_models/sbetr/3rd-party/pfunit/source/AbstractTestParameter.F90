module AbstractTestParameter_mod
   implicit none
   private

   public :: AbstractTestParameter

   type, abstract :: AbstractTestParameter
   contains
      procedure(toString), deferred :: toString
      procedure :: toStringActual
   end type AbstractTestParameter

   abstract interface
      function toString(this) result(string)
         import AbstractTestParameter
         class (AbstractTestParameter), intent(in) :: this
         character(:), allocatable :: string
      end function toString
   end interface

contains

   function toStringActual(this) result(string)
      class (AbstractTestParameter), intent(in) :: this
      character(:), allocatable :: string
      string = this%toString()
   end function toStringActual

end module AbstractTestParameter_mod
