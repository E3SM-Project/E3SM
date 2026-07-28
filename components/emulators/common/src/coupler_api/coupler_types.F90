module coupler_types
   use iso_c_binding
   implicit none

   type, bind(c) :: emulator_create_cfg
      integer(c_int) :: f_comm
      integer(c_int) :: comp_id
      integer(c_int) :: run_type
      integer(c_int) :: start_ymd
      integer(c_int) :: start_tod
      integer(c_int) :: case_start_ymd
      integer(c_int) :: case_start_tod
      type(c_ptr)    :: input_file  ! C char*
      type(c_ptr)    :: log_file    ! C char*
      type(c_ptr)    :: calendar    ! C char*
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

   type, bind(c) :: coupling_desc
      type(c_ptr)    :: import_data  ! double*
      type(c_ptr)    :: export_data  ! double*
      integer(c_int) :: num_imports
      integer(c_int) :: num_exports
      integer(c_int) :: field_size
   end type coupling_desc

   type, bind(c) :: field_attributes
      ! these are all C char *
      type(c_ptr) name; 
      type(c_ptr) long_name; 
      type(c_ptr) standard_name; 
      type(c_ptr) units; 
   end type

   type, bind(c) :: registered_field_desc
      type(c_ptr) role; 
      type(c_ptr) component; 
      type(c_ptr) attributes; 
      integer(c_size_t) size; 
      type(c_ptr) data; ! double*
   end type

   public :: create_config, create_grid_desc, create_coupler_desc
   public :: coupling_desc, emulator_grid_desc, emulator_create_cfg
   public :: field_attributes, registered_field_desc

contains

   subroutine to_c_string(str, cstr)
      character(len=*), intent(in) :: str
      character(kind=c_char), allocatable, target, intent(out) :: cstr(:)

      integer :: i, n

      n = len_trim(str)
      allocate (cstr(n + 1))

      do i = 1, n
         cstr(i) = str(i:i)
      end do
      cstr(n + 1) = c_null_char
   end subroutine to_c_string

   type(emulator_create_cfg) &
      function create_config(f_comm, comp_id, run_type, &
                             start_ymd, start_tod, case_start_ymd, case_start_tod, &
                             input_file, log_file, calendar) result(cfg)
      integer(c_int), intent(in) :: f_comm
      integer(c_int), intent(in) :: comp_id
      integer(c_int), intent(in) :: run_type
      integer(c_int), intent(in) :: start_ymd
      integer(c_int), intent(in) :: start_tod
      integer(c_int), intent(in), optional :: case_start_ymd
      integer(c_int), intent(in), optional :: case_start_tod
      character(kind=c_char), intent(in), target   :: input_file(*)
      character(kind=c_char), intent(in), target   :: log_file(*)
      character(kind=c_char), intent(in), target   :: calendar(*)

      cfg%f_comm = f_comm
      cfg%comp_id = comp_id
      cfg%run_type = run_type
      cfg%start_ymd = start_ymd
      cfg%start_tod = start_tod
      if (present(case_start_ymd)) then
         cfg%case_start_ymd = case_start_ymd
      else
         cfg%case_start_ymd = -1
      end if
      if (present(case_start_tod)) then
         cfg%case_start_tod = case_start_tod
      else
         cfg%case_start_tod = -1
      end if
      cfg%input_file = c_loc(input_file(1))
      cfg%log_file = c_loc(log_file(1))
      cfg%calendar = c_loc(calendar(1))

   end function create_config

   type(emulator_grid_desc) &
      function create_grid_desc( &
      grid_type, nx, ny, num_local_cols, num_global_cols, &
      col_gids, lat, lon, area) result(grid)
      integer(c_int), intent(in) :: grid_type
      integer(c_int), intent(in) :: nx
      integer(c_int), intent(in) :: ny
      integer(c_int), intent(in) :: num_local_cols
      integer(c_int), intent(in) :: num_global_cols
      integer(c_int), intent(in), pointer :: col_gids(:)
      real(c_double), intent(in), pointer :: lat(:)
      real(c_double), intent(in), pointer :: lon(:)
      real(c_double), intent(in), pointer :: area(:)

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

   type(coupling_desc) &
      function create_coupler_desc(import_data, export_data, num_imports, &
                                   num_exports, field_size) result(cpl)
      real(c_double), INTENT(IN), pointer :: import_data(:)
      real(c_double), INTENT(IN), pointer :: export_data(:)
      integer(c_int) :: num_imports
      integer(c_int) :: num_exports
      integer(c_int) :: field_size
      cpl%import_data = c_loc(import_data)
      cpl%export_data = c_loc(export_data)
      cpl%num_imports = num_imports
      cpl%num_exports = num_exports
      cpl%field_size = field_size
   end function create_coupler_desc

   type(field_attributes) &
      function make_field_attributes(name, long_name, standard_name, units) result(attr)
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: long_name
      character(len=*), intent(in) :: standard_name
      character(len=*), intent(in) :: units

      character(kind=c_char), allocatable, target :: name_c(:)
      character(kind=c_char), allocatable, target :: long_name_c(:)
      character(kind=c_char), allocatable, target :: standard_name_c(:)
      character(kind=c_char), allocatable, target :: units_c(:)

      call to_c_string(name, name_c)
      call to_c_string(long_name, long_name_c)
      call to_c_string(standard_name, standard_name_c)
      call to_c_string(units, units_c)

      attr%name = c_loc(name_c(1))
      attr%long_name = c_loc(long_name_c(1))
      attr%standard_name = c_loc(standard_name_c(1))
      attr%units = c_loc(units_c(1))

   end function make_field_attributes

   type(registered_field_desc) &
      function make_registered_field_desc(role, component, attributes, size, data) result(desc)
      character(len=*), intent(in) :: role
      character(len=*), intent(in) :: component
      character(len=*), intent(in) :: attributes
      integer, intent(in) :: size
      real(c_double), pointer, intent(in) :: data(:) ! double*
      !! local c versions of strings
      character(kind=c_char), allocatable, target :: role_c(:)
      character(kind=c_char), allocatable, target :: component_c(:)
      character(kind=c_char), allocatable, target :: attributes_c(:)

      call to_c_string(role, role_c)
      call to_c_string(component, component_c)
      call to_c_string(attributes, attributes_c)

      desc%role = c_loc(role_c(1))
      desc%component = c_loc(component_c(1))
      desc%attributes = c_loc(attributes_c(1))

      desc%size = int(size, c_size_t); 
      desc%data = c_loc(data(1)); 
   end function make_registered_field_desc

end module coupler_types
