program test_emulator_factory
  use iso_c_binding
  use emulator_f_api
  use emulator_f2c_api
  use emulator_handle_mod, only : emulator_handle
  use mpi

  implicit none

  type(emulator_create_cfg)    :: cfg
  type(emulator_grid_desc)     :: grid
  type(emulator_coupling_desc) :: cpl
  type(c_ptr)                  :: handle
  integer(c_int)               :: dt
  integer                      :: ierr, nprocs
  integer(c_int)               :: fcomm
  integer(c_int), parameter    :: nx = 1_c_int, ny = 1_c_int
  integer(c_int), parameter    :: ncols = 1_c_int
  integer(c_int), parameter :: num_global_cols = 1_c_int

  integer(c_int) :: num_local_cols, field_size, num_imports, num_exports
  integer(c_int), allocatable, target:: col_gids(:)
  real(c_double),allocatable, target:: lat(:), lon(:), area(:)
  real(c_double), allocatable, target:: import_buf(:), export_buf(:)
  real(c_double), pointer :: lat_ptr(:),lon_ptr(:), imp_ptr(:), exp_ptr(:),area_ptr(:)
  integer(c_int), pointer :: gids_ptr(:)

  !----------------------------------------
  ! MPI init (if needed)
  !----------------------------------------
  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) then
     print *, "MPI_Init failed"
     stop 1
  end if
  fcomm = MPI_COMM_WORLD
  call mpi_comm_size(fcomm, nprocs, ierr)

   num_local_cols = num_global_cols/nprocs
   num_imports = 1_c_int
   num_exports = 1_c_int
   field_size = 1_c_int

   allocate(lat(num_local_cols), lon(num_local_cols), col_gids(num_local_cols))
   allocate(area(num_local_cols),import_buf(num_imports*field_size), export_buf(num_exports*field_size))
  !----------------------------------------
  ! create config
  !----------------------------------------
  cfg = create_config(f_comm=fcomm,comp_id=1_c_int,run_type=0_c_int,&
            start_ymd=20000101_c_int, start_tod=0_c_int,&
            input_file='test'//c_null_char, log_file="test_log"//c_null_char)
  
  block
   integer(c_int) :: grid_type = 0_c_int, i
   lat = [(real(i,kind=c_double), i=1,num_local_cols)]
   lon = lat
   area = lat
   lat_ptr => lat; lon_ptr => lon; area_ptr => area

   col_gids = [(i, i=1,num_global_cols) ]
   gids_ptr => col_gids

   grid = create_grid_desc(grid_type=grid_type,nx=nx,ny=ny,&
         num_local_cols=num_local_cols, num_global_cols=num_global_cols,&
         col_gids=gids_ptr, lat=lat_ptr, lon=lon_ptr, area=area_ptr)

  import_buf(:) = [(real(i*i,kind=c_double), i=1,num_imports*field_size)]
  export_buf = [(real(i*i,kind=c_double), i=1,num_exports*field_size)]
  end block
  imp_ptr => import_buf; exp_ptr => export_buf
  cpl = create_coupler_desc(import_data=imp_ptr,export_data=exp_ptr,&
               num_imports=num_imports,num_exports=num_exports,field_size=field_size)
  !----------------------------------------
  ! Call emulator_create("atm")
  !----------------------------------------
  handle = emulator_create("atm"//c_null_char, cfg)

  if (.not. c_associated(handle)) then
     print *, "ERROR: emulator_create returned NULL handle for 'atm'"
     stop 1
  else
     print *, "OK: emulator_create returned non-null handle for 'atm'"
  end if

  call emulator_set_grid_data(handle, grid)
  call emulator_setup_coupling(handle, cpl)

  !----------------------------------------
  ! Init / run / finalize
  !----------------------------------------
  call emulator_init(handle)

  dt = 3600_c_int
  call emulator_run(handle, dt)
  call emulator_finalize(handle)

  print *, "OK: basic emulator lifecycle completed"

#ifdef USE_MPI
  call MPI_Finalize(ierr)
#endif

end program test_emulator_factory
