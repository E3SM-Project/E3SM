module emulator_f2c_api
   use, intrinsic :: iso_c_binding
   use emulator_f_api, only: emulator_coupling_desc, emulator_grid_desc, emulator_create_cfg
   implicit none

   public
   interface
      function emulator_create(kind, cfg) result(handle) bind(c)
         import :: c_ptr, c_char, emulator_create_cfg
         character(kind=c_char), intent(in) :: kind(*)
         type(emulator_create_cfg), intent(in) :: cfg
         type(c_ptr) :: handle
      end function emulator_create

      subroutine emulator_set_grid_data(handle, grid) bind(c)
         import :: c_ptr, emulator_grid_desc
         type(c_ptr), value, intent(in) :: handle
         type(emulator_grid_desc), intent(in) :: grid
      end subroutine emulator_set_grid_data

      subroutine emulator_setup_coupling(handle, cpl) bind(c)
         import :: c_ptr, emulator_coupling_desc
         type(c_ptr), value, intent(in) :: handle
         type(emulator_coupling_desc), intent(in) :: cpl
      end subroutine emulator_setup_coupling

      subroutine emulator_init(handle) bind(c)
         import :: c_ptr
         type(c_ptr), value, intent(in) :: handle
      end subroutine emulator_init

      subroutine emulator_run(handle, dt) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: dt
      end subroutine emulator_run

      subroutine emulator_finalize(handle) bind(c)
         import :: c_ptr
         type(c_ptr), value, intent(in) :: handle
      end subroutine emulator_finalize

      subroutine emulator_print_info(handle) bind(c)
         import :: c_ptr
         type(c_ptr), value, intent(in) :: handle
      end subroutine emulator_print_info

      subroutine emulator_init_coupling_indices(handle,import_fields, export_fields) bind(c)
         import :: c_char, c_ptr
         type(c_ptr), value, intent(in) :: handle
         type(c_ptr), value :: import_fields ! char*
         type(c_ptr), value :: export_fields ! char*
      end subroutine emulator_init_coupling_indices

    function emulator_get_num_local_cols(handle) result(n) bind(c)
      import :: c_ptr, c_int
      type(c_ptr),              value :: handle
      integer(c_int)                 :: n
    end function emulator_get_num_local_cols

    function emulator_get_num_global_cols(handle) result(n) bind(c)
      import :: c_ptr, c_int
      type(c_ptr),              value :: handle
      integer(c_int)                 :: n
    end function emulator_get_num_global_cols

    function emulator_get_nx(handle) result(nx) bind(c)
      import :: c_ptr, c_int
      type(c_ptr),              value :: handle
      integer(c_int)                 :: nx
    end function emulator_get_nx

    function emulator_get_ny(handle) result(ny) bind(c)
      import :: c_ptr, c_int
      type(c_ptr),              value :: handle
      integer(c_int)                 :: ny
    end function emulator_get_ny

    subroutine emulator_get_local_col_gids(handle, gids) bind(c)
      import :: c_ptr, c_int
      type(c_ptr),              value :: handle
      integer(c_int), dimension(*), intent(out) :: gids
    end subroutine emulator_get_local_col_gids

    subroutine emulator_get_cols_latlon(handle, lat, lon) bind(c)
      import :: c_ptr, c_double
      type(c_ptr),              value :: handle
      real(c_double), dimension(*), intent(out) :: lat
      real(c_double), dimension(*), intent(out) :: lon
    end subroutine emulator_get_cols_latlon

    subroutine emulator_get_cols_area(handle, area) bind(c)
      import :: c_ptr, c_double
      type(c_ptr),              value :: handle
      real(c_double), dimension(*), intent(out) :: area
    end subroutine emulator_get_cols_area
   end interface

end module emulator_f2c_api
