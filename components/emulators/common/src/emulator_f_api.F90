module emulator_f_api
  use iso_c_binding
  implicit none

  type, bind(c) :: emulator_create_cfg
     integer(c_int) :: f_comm
     integer(c_int) :: comp_id
     integer(c_int) :: run_type
     integer(c_int) :: start_ymd
     integer(c_int) :: start_tod
     type(c_ptr)    :: input_file  ! C char*
     type(c_ptr)    :: log_file    ! C char*
  end type emulator_create_cfg

  type, bind(c) :: emulator_grid_desc
     integer(c_int) :: grid_type
     integer(c_int) :: nx
     integer(c_int) :: ny
     integer(c_int) :: num_local_cols
     integer(c_int) :: num_global_cols
     type(c_ptr)    :: col_gids  ! int*
     type(c_ptr)    :: lat       ! double*
     type(c_ptr)    :: lon       ! double*
     type(c_ptr)    :: area      ! double*
  end type emulator_grid_desc

  type, bind(c) :: emulator_coupling_desc
     type(c_ptr)    :: import_data  ! double*
     type(c_ptr)    :: export_data  ! double*
     integer(c_int) :: num_imports
     integer(c_int) :: num_exports
     integer(c_int) :: field_size
  end type emulator_coupling_desc

  public :: create_config, create_grid_desc, create_coupler_desc

contains

   type(emulator_create_cfg) &
   function create_config(f_comm, comp_id, run_type,&
            start_ymd, start_tod, &
            input_file, log_file) result(cfg)
      integer(c_int), intent(in) :: f_comm
      integer(c_int), intent(in) :: comp_id
      integer(c_int), intent(in) :: run_type
      integer(c_int), intent(in) :: start_ymd
      integer(c_int), intent(in) :: start_tod
      character(kind=c_char), intent(in), target   :: input_file(*)
      character(kind=c_char), intent(in), target   :: log_file(*)

      cfg%f_comm=f_comm
      cfg%comp_id=comp_id
      cfg%run_type=run_type
      cfg%start_ymd=start_ymd
      cfg%start_tod=start_tod
      cfg%input_file=c_loc(input_file(1))
      cfg%log_file=c_loc(log_file(1))

   end function create_config

   type(emulator_grid_desc)&
   function create_grid_desc( &
      grid_type, nx, ny, num_local_cols, num_global_cols,&
      col_gids, lat, lon, area) result(grid)
     integer(c_int),intent(in) :: grid_type
     integer(c_int),intent(in) :: nx
     integer(c_int),intent(in) :: ny
     integer(c_int),intent(in) :: num_local_cols
     integer(c_int),intent(in) :: num_global_cols
     integer(c_int),intent(in), pointer :: col_gids(:)
     real(c_double),intent(in), pointer :: lat(:)
     real(c_double),intent(in), pointer :: lon(:)
     real(c_double),intent(in), pointer :: area(:)

      grid%grid_type = grid_type
      grid%nx = nx
      grid%ny = ny
      grid%num_local_cols = num_local_cols
      grid%num_global_cols = num_global_cols
      grid%col_gids = c_loc(col_gids)
      grid%lat = c_loc(lat)
      grid%lon = c_loc(lon)
      grid%area = c_loc(area)
   end function create_grid_desc

   type(emulator_coupling_desc)&
      function create_coupler_desc(import_data,export_data,num_imports,&
         num_exports,field_size) result(cpl)
        real(c_double),INTENT(IN), pointer :: import_data(:)
        real(c_double),INTENT(IN), pointer :: export_data(:)
        integer(c_int) :: num_imports
        integer(c_int) :: num_exports
        integer(c_int) :: field_size
        cpl%import_data = c_loc(import_data)
        cpl%export_data = c_loc(export_data)
        cpl%num_imports = num_imports
        cpl%num_exports = num_exports
        cpl%field_size = field_size
     end function create_coupler_desc

end module emulator_f_api
