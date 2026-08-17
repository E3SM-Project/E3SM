module omega_f2cxx_mod

   implicit none
   public

   ! All type(c_ptr) dummy args below use VALUE: c_loc() already yields the
   ! address to pass, so VALUE hands that address straight through as a
   ! plain C/C++ pointer, matching the raw-pointer C++ signatures. None of
   ! these callees reseat the pointer, so a reference-to-pointer is not
   ! needed on the C++ side.

   interface
      subroutine omega_ocn_init1( &
         f_comm, &
         ocn_id, &
         yaml_config_name, &
         ocn_log_name, &
         start_type, &
         calendar_name, &
         run_start_ymd, &
         run_start_tod, &
         coupler_time_step, &
         n_coupler_imports, &
         n_coupler_exports, &
         n_omega_imports, &
         n_omega_exports, &
         import_field_names, &
         export_field_names, &
         import_field_indices, &
         export_field_indices, &
         cpl_x2o_field_names, &
         cpl_o2x_field_names) bind(c)

         use, intrinsic :: iso_c_binding, only: c_int, c_char, c_ptr

         implicit none

         integer(kind=c_int), value, intent(in) :: &
            f_comm, &
            ocn_id, &
            start_type, &
            run_start_ymd, &
            run_start_tod, &
            coupler_time_step, &
            n_coupler_imports, &
            n_coupler_exports, &
            n_omega_imports, &
            n_omega_exports

         character(kind=c_char), target, intent(in) :: &
            yaml_config_name, ocn_log_name, calendar_name

         type(c_ptr), value, intent(in) :: &
            import_field_names, &
            export_field_names, &
            import_field_indices, &
            export_field_indices, &
            cpl_x2o_field_names, &
            cpl_o2x_field_names

      end subroutine omega_ocn_init1

      subroutine omega_ocn_init2(cpl_to_ocn_data, ocn_to_cpl_data &
                                 ) bind(c)

         use, intrinsic :: iso_c_binding, only: c_ptr

         implicit none

         type(c_ptr), value, intent(in) :: &
            cpl_to_ocn_data, ocn_to_cpl_data

      end subroutine omega_ocn_init2

      subroutine omega_ocn_run(write_restart) bind(c)

         use, intrinsic :: iso_c_binding, only: c_bool

         implicit none

         logical(kind=c_bool), value, intent(in) :: write_restart

      end subroutine omega_ocn_run

      subroutine omega_ocn_finalize() bind(c)

         implicit none

      end subroutine omega_ocn_finalize

      function omega_get_moab_pid() result(moab_pid) bind(c)

         use, intrinsic :: iso_c_binding, only: c_int

         implicit none

         integer(kind=c_int) :: moab_pid

      end function omega_get_moab_pid

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
