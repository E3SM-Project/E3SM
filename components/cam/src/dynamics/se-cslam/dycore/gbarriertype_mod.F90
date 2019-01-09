module gbarriertype_mod
  use ISO_C_Binding, only: C_ptr

   type, public :: gbarrier_t
     type (C_ptr) :: c_barrier
   end type gbarrier_t

end module gbarriertype_mod
