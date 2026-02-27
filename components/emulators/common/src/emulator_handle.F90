module emulator_handle_mod
  use, intrinsic :: iso_c_binding
  use emulator_f2c_api
  implicit none

  type :: emulator_handle
     type(c_ptr) :: h = c_null_ptr
   contains
     procedure :: initialize
     procedure :: run
     !! setters
     procedure :: set_grid_data
     procedure :: setup_coupling

     !! getters
     procedure :: num_local_cols  => emulator_num_local_cols
     procedure :: num_global_cols => emulator_num_global_cols
     procedure :: nx              => emulator_nx
     procedure :: ny              => emulator_ny
     procedure :: get_local_col_gids
     procedure :: get_cols_latlon
     procedure :: get_cols_area
     !! diags
     procedure :: print_info
     !! Clean-up
     final  :: destroy_emulator_handle
     procedure :: finalize => emulator_finalize_wrapper
  end type emulator_handle

contains

  subroutine set_grid_data(self, grid)
    class(emulator_handle),      intent(inout) :: self
    type(emulator_grid_desc),    intent(in)    :: grid
    call emulator_set_grid_data(self%h, grid)
  end subroutine set_grid_data

  subroutine setup_coupling(self, cpl)
    class(emulator_handle),         intent(inout) :: self
    type(emulator_coupling_desc),   intent(in)    :: cpl
    call emulator_setup_coupling(self%h, cpl)
  end subroutine setup_coupling

  subroutine initialize(self)
    class(emulator_handle), intent(inout) :: self
    call emulator_init(self%h)
  end subroutine initialize

  subroutine run(self, dt)
    class(emulator_handle), intent(inout) :: self
    integer(c_int),         intent(in)    :: dt
    call emulator_run(self%h, dt)
  end subroutine run

  subroutine destroy_emulator_handle(self)
    type(emulator_handle), intent(inout) :: self
    call emulator_finalize_wrapper(self)
  end subroutine destroy_emulator_handle

  subroutine emulator_finalize_wrapper(self)
    class(emulator_handle), intent(inout) :: self
    if (c_associated(self%h)) then
      call emulator_finalize(self%h)
      self%h = c_null_ptr
    end if
  end subroutine emulator_finalize_wrapper

  integer function emulator_num_local_cols(self) result(n)
    class(emulator_handle), intent(in) :: self
    n = emulator_get_num_local_cols(self%h)
  end function emulator_num_local_cols

  integer function emulator_num_global_cols(self) result(n)
    class(emulator_handle), intent(in) :: self
    n = emulator_get_num_global_cols(self%h)
  end function emulator_num_global_cols

  integer function emulator_nx(self) result(nx_val)
    class(emulator_handle), intent(in) :: self
    nx_val = emulator_get_nx(self%h)
  end function emulator_nx

  integer function emulator_ny(self) result(ny_val)
    class(emulator_handle), intent(in) :: self
    ny_val = emulator_get_ny(self%h)
  end function emulator_ny

  subroutine get_local_col_gids(self, gids)
    class(emulator_handle), intent(in) :: self
    integer,               intent(out) :: gids(:)
    call emulator_get_local_col_gids(self%h, gids)
  end subroutine get_local_col_gids

  subroutine get_cols_latlon(self, lat, lon)
    class(emulator_handle), intent(in) :: self
    real(c_double),              intent(out) :: lat(:), lon(:)
    call emulator_get_cols_latlon(self%h, lat, lon)
  end subroutine get_cols_latlon

  subroutine get_cols_area(self, area)
    class(emulator_handle), intent(in) :: self
    real(c_double), intent(out) :: area(:)
    call emulator_get_cols_area(self%h, area)
  end subroutine get_cols_area

  subroutine print_info(self)
    class(emulator_handle), intent(in) :: self
    call emulator_print_info(self%h)
  end subroutine print_info

  subroutine init_coupling_indices(self, import_fields,export_fields)
    class(emulator_handle), intent(in) :: self
    character(kind=c_char), intent(in), target :: import_fields(*)
    character(kind=c_char), intent(in), target :: export_fields(*)
    call emulator_init_coupling_indices(self%h,&
          c_loc(import_fields(1)), &
          c_loc(export_fields(1)))
  end subroutine init_coupling_indices

end module emulator_handle_mod
