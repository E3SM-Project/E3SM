module scream_rrtmgp_interface_mod
   use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_double, c_bool, C_NULL_CHAR, c_float
   implicit none
   public rrtmgp_init_f90
   public rrtmgp_main_f90
   public rrtmgp_finalize_f90

contains

   subroutine rrtmgp_init_f90()
   end subroutine rrtmgp_init_f90

   subroutine rrtmgp_main_f90()
   end subroutine rrtmgp_main_f90

   subroutine rrtmgp_finalize_f90()
   end subroutine rrtmgp_finalize_f90

end module scream_rrtmgp_interface_mod
