module omega_f2cxx_mod

   implicit none
   public

   interface
      subroutine omega_ocn_init( &
         f_comm, &
         ocn_id, &
         yaml_config_name, &
         ocn_log_name, &
         calendar_name, &
         run_start_ymd, &
         run_start_tod) bind(c)

         use, intrinsic :: iso_c_binding, only: c_int, c_char

         implicit none

         integer(kind=c_int), value, intent(in) :: &
            f_comm, ocn_id, run_start_ymd, run_start_tod
         character(kind=c_char), target, intent(in) :: &
            yaml_config_name, ocn_log_name, calendar_name

      end subroutine omega_ocn_init
   end interface
end module omega_f2cxx_mod
