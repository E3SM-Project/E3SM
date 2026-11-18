module harness_handle_mod
   use iso_c_binding, only: c_ptr
   implicit none
   private

   type, public :: harness_handle
      type(c_ptr) :: handle = c_ptr
   contains
      procedure  :: init
   end type harness_handle

contains

end module harness_handle_mod
