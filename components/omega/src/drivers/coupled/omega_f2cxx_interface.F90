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

      function omega_get_ncells_local() result(ncells_local) bind(c)

         use, intrinsic :: iso_c_binding, only: c_int

         implicit none

         integer(kind=c_int) :: ncells_local

      end function omega_get_ncells_local

      function omega_get_ncells_global() result(ncells_global) bind(c)

         use, intrinsic :: iso_c_binding, only: c_int

         implicit none

         integer(kind=c_int) :: ncells_global

      end function omega_get_ncells_global

      subroutine omega_get_index_to_cell_id( &
         local_index_to_cell_id &
         ) bind(c)

         use, intrinsic :: iso_c_binding, only: c_int

         implicit none

         ! allow(assumed-size)
         integer(kind=c_int), intent(out) :: local_index_to_cell_id(*)

      end subroutine omega_get_index_to_cell_id

      subroutine omega_get_lonlat_cell(lon_cell, lat_cell) bind(c)

         use, intrinsic :: iso_c_binding, only: c_double

         implicit none

         ! allow(assumed-size)
         real(kind=c_double), intent(out) :: lon_cell(*), lat_cell(*)
      end subroutine omega_get_lonlat_cell

      subroutine omega_get_area_cell(area_cell) bind(c)

         use, intrinsic :: iso_c_binding, only: c_double

         implicit none

         ! allow(assumed-size)
         real(kind=c_double), intent(out) :: area_cell(*)
      end subroutine omega_get_area_cell
   end interface
end module omega_f2cxx_mod
