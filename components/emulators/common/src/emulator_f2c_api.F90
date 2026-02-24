module emulator_f2c_api
  use, intrinsic :: iso_c_binding
  use emulator_f_api, only: emulator_coupling_desc, emulator_grid_desc, emulator_create_cfg
  implicit none
  private

  public :: emulator_create, emulator_set_grid_data, emulator_setup_coupling
  public :: emulator_init, emulator_run, emulator_finalize

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
       use emulator_f_api, only: emulator_coupling_desc
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
  end interface

end module emulator_f2c_api
