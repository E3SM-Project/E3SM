module emulator_atm_f2c
   !--------------------------------------------------------------------------
   ! Fortran interface declarations for the C++ EmulatorAtm class.
   !
   ! Provides iso_c_binding interfaces that allow the Fortran MCT wrapper
   ! (atm_comp_mct.F90) to call the C++ EmulatorAtm implementation.
   !--------------------------------------------------------------------------

   use iso_c_binding
   implicit none
   private

   public :: emulator_atm_create_instance
   public :: emulator_atm_set_grid_data
   public :: emulator_atm_init_coupling_indices
   public :: emulator_atm_setup_coupling
   public :: emulator_atm_init
   public :: emulator_atm_run
   public :: emulator_atm_finalize
   public :: emulator_atm_get_num_local_cols
   public :: emulator_atm_get_num_global_cols
   public :: emulator_atm_get_nx
   public :: emulator_atm_get_ny
   public :: emulator_atm_get_local_cols_gids
   public :: emulator_atm_get_cols_latlon
   public :: emulator_atm_get_cols_area

   interface

      subroutine emulator_atm_create_instance(f_comm, comp_id, &
         input_file, log_file, run_type, start_ymd, start_tod) &
         bind(c)
         use iso_c_binding
         integer(c_int), value, intent(in) :: f_comm
         integer(c_int), value, intent(in) :: comp_id
         character(kind=c_char), intent(in) :: input_file(*)
         character(kind=c_char), intent(in) :: log_file(*)
         integer(c_int), value, intent(in) :: run_type
         integer(c_int), value, intent(in) :: start_ymd
         integer(c_int), value, intent(in) :: start_tod
      end subroutine emulator_atm_create_instance

      subroutine emulator_atm_set_grid_data(nx, ny, &
         num_local_cols, num_global_cols, &
         col_gids, lat, lon, area) bind(c)
         use iso_c_binding
         integer(c_int), value, intent(in) :: nx
         integer(c_int), value, intent(in) :: ny
         integer(c_int), value, intent(in) :: num_local_cols
         integer(c_int), value, intent(in) :: num_global_cols
         type(c_ptr), value, intent(in) :: col_gids
         type(c_ptr), value, intent(in) :: lat
         type(c_ptr), value, intent(in) :: lon
         type(c_ptr), value, intent(in) :: area
      end subroutine emulator_atm_set_grid_data

      subroutine emulator_atm_init_coupling_indices( &
         export_fields, import_fields) bind(c)
         use iso_c_binding
         character(kind=c_char), intent(in) :: export_fields(*)
         character(kind=c_char), intent(in) :: import_fields(*)
      end subroutine emulator_atm_init_coupling_indices

      subroutine emulator_atm_setup_coupling(import_data, &
         export_data, num_imports, num_exports, field_size) &
         bind(c)
         use iso_c_binding
         type(c_ptr), value, intent(in) :: import_data
         type(c_ptr), value, intent(in) :: export_data
         integer(c_int), value, intent(in) :: num_imports
         integer(c_int), value, intent(in) :: num_exports
         integer(c_int), value, intent(in) :: field_size
      end subroutine emulator_atm_setup_coupling

      subroutine emulator_atm_init() bind(c)
      end subroutine emulator_atm_init

      subroutine emulator_atm_run(dt) bind(c)
         use iso_c_binding
         integer(c_int), value, intent(in) :: dt
      end subroutine emulator_atm_run

      subroutine emulator_atm_finalize() bind(c)
      end subroutine emulator_atm_finalize

      function emulator_atm_get_num_local_cols() result(n) &
         bind(c)
         use iso_c_binding
         integer(c_int) :: n
      end function emulator_atm_get_num_local_cols

      function emulator_atm_get_num_global_cols() result(n) &
         bind(c)
         use iso_c_binding
         integer(c_int) :: n
      end function emulator_atm_get_num_global_cols

      function emulator_atm_get_nx() result(n) bind(c)
         use iso_c_binding
         integer(c_int) :: n
      end function emulator_atm_get_nx

      function emulator_atm_get_ny() result(n) bind(c)
         use iso_c_binding
         integer(c_int) :: n
      end function emulator_atm_get_ny

      subroutine emulator_atm_get_local_cols_gids(gids) bind(c)
         use iso_c_binding
         type(c_ptr), value, intent(in) :: gids
      end subroutine emulator_atm_get_local_cols_gids

      subroutine emulator_atm_get_cols_latlon(lat, lon) bind(c)
         use iso_c_binding
         type(c_ptr), value, intent(in) :: lat
         type(c_ptr), value, intent(in) :: lon
      end subroutine emulator_atm_get_cols_latlon

      subroutine emulator_atm_get_cols_area(area) bind(c)
         use iso_c_binding
         type(c_ptr), value, intent(in) :: area
      end subroutine emulator_atm_get_cols_area

   end interface

end module emulator_atm_f2c
