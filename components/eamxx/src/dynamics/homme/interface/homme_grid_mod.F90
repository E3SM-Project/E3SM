module homme_grid_mod
  use iso_c_binding,     only: c_ptr, c_f_pointer, c_int, c_bool, c_double
  use parallel_mod,      only: abortmp
  use kinds,             only: iulog
  use homme_context_mod, only: par, masterproc

  implicit none
  private

  ! Routines that modify state
  public :: init_grids_f90
  public :: finalize_geometry_f90

  ! Routines to get information
  ! public :: get_cols_gids_f90, get_cols_indices_f90
  public :: get_num_local_columns_f90, get_num_global_columns_f90
  public :: get_num_local_elems_f90, get_num_global_elems_f90
  public :: get_np_f90, get_nlev_f90

contains

  subroutine check_grids_inited (expected)
    use parallel_mod,      only: abortmp
    use homme_context_mod, only: is_geometry_inited
    !
    ! Input(s)
    !
    logical, intent(in) :: expected

    if (expected) then
      if (.not. is_geometry_inited) then
        call abortmp ("Error! Geometry has not been initialized yet (or has already been finalized).\n")
      endif
    else
      if (is_geometry_inited) then
        call abortmp ("Error! Geometry was already initialized.\n")
      endif
    endif
  end subroutine check_grids_inited

  subroutine init_grids_f90 (pg_types_ptr,num_pg_types) bind(c)
    use dyn_grid_mod,      only: dyn_grid_init
    use phys_grid_mod,     only: phys_grids_init
    use homme_context_mod, only: is_geometry_inited
    !
    ! Input(s)
    !
    integer (kind=c_int), intent(in), value :: num_pg_types
    type(c_ptr), intent(in) :: pg_types_ptr
    !
    ! Local(s)
    !
    integer (kind=c_int), pointer :: pg_types(:)

    call c_f_pointer(pg_types_ptr, pg_types, [num_pg_types])

    ! We don't expect this to be called twice
    call check_grids_inited(.false.)

    ! Init SE grid struff in Homme
    call dyn_grid_init ()

    ! Init physics grids
    call phys_grids_init (pg_types)

    is_geometry_inited = .true.
  end subroutine init_grids_f90

  subroutine cleanup_grid_init_data_f90 () bind(c)
    use phys_grid_mod, only: pg_cleanup_grid_init_data=>cleanup_grid_init_data
    use dyn_grid_mod,  only: dg_cleanup_grid_init_data=>cleanup_grid_init_data

    ! Do not call this if grids were not inited!
    call check_grids_inited(.true.)

    ! Clean up in dyn and phys grid mods
    call dg_cleanup_grid_init_data()
    call pg_cleanup_grid_init_data ()
  end subroutine cleanup_grid_init_data_f90

  subroutine finalize_geometry_f90 () bind(c)
    use homme_context_mod, only: is_geometry_inited
    use phys_grid_mod,     only: finalize_phys_grid

    ! Don't finalize what you didn't initialize.
    call check_grids_inited(.true.)

    call finalize_phys_grid ()

    is_geometry_inited = .false.
  end subroutine finalize_geometry_f90

  ! Note: dg_grid=.true. is to request dofs for a Discontinuous Galerkin grid,
  !       that is, corresponding edge dofs on bordering elems have different gids.
  !       If dg_grid=.false., the shared dofs have the same gid.
  subroutine get_dyn_grid_data_f90 (dg_gids_ptr, cg_gids_ptr, elgpgp_ptr, lat_ptr, lon_ptr) bind(c)
    use dimensions_mod,    only: nelemd, np
    use dyn_grid_mod,      only: get_my_dyn_data
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: dg_gids_ptr, cg_gids_ptr, elgpgp_ptr, lat_ptr, lon_ptr
    !
    ! Local(s)
    !
    real(kind=c_double), pointer :: lat (:,:,:), lon(:,:,:)
    integer(kind=c_int), pointer :: cg_gids (:), dg_gids(:), elgpgp(:,:)

    ! Sanity check
    call check_grids_inited(.true.)

    call c_f_pointer (dg_gids_ptr, dg_gids, [nelemd*np*np])
    call c_f_pointer (cg_gids_ptr, cg_gids, [nelemd*np*np])
    call c_f_pointer (elgpgp_ptr,  elgpgp,  [3,nelemd*np*np])
    call c_f_pointer (lat_ptr,     lat,     [np,np,nelemd])
    call c_f_pointer (lon_ptr,     lon,     [np,np,nelemd])

    call get_my_dyn_data (dg_gids, cg_gids, elgpgp, lat, lon)
  end subroutine get_dyn_grid_data_f90

  subroutine get_phys_grid_data_f90 (pg_type, gids_ptr, lat_ptr, lon_ptr, area_ptr) bind(c)
    use phys_grid_mod, only: get_my_phys_data
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: gids_ptr, lat_ptr, lon_ptr, area_ptr
    integer (kind=c_int), intent(in) :: pg_type
    !
    ! Local(s)
    !
    integer :: ncols
    real(kind=c_double), pointer :: lat(:), lon(:), area(:)
    integer(kind=c_int), pointer :: gids(:)

    ! Sanity check
    call check_grids_inited(.true.)

    ncols = get_num_local_columns_f90(mod(pg_type,10))

    call c_f_pointer (gids_ptr, gids, [ncols])
    call c_f_pointer (lat_ptr,  lat,  [ncols])
    call c_f_pointer (lon_ptr,  lon,  [ncols])
    call c_f_pointer (area_ptr, area, [ncols])

    ! Retrieve all data
    call get_my_phys_data (gids, lat, lon, area, pg_type)

  end subroutine get_phys_grid_data_f90

  function get_num_local_columns_f90 (pg_type) result (ncols) bind(c)
    use phys_grid_mod,     only: get_num_local_columns
    !
    ! Input(s)
    !
    integer (kind=c_int), intent(in), value :: pg_type
    !
    ! Local(s)
    !
    integer (kind=c_int) :: ncols

    ! Sanity check
    call check_grids_inited(.true.)

    ncols = get_num_local_columns(mod(pg_type,10))
  end function get_num_local_columns_f90

  function get_num_global_columns_f90 (pg_type) result (ncols) bind(c)
    use phys_grid_mod, only: get_num_global_columns
    !
    ! Input(s)
    !
    integer (kind=c_int), intent(in), value :: pg_type
    !
    ! Local(s)
    !
    integer (kind=c_int) :: ncols

    ! Sanity check
    call check_grids_inited(.true.)

    ncols = get_num_global_columns(mod(pg_type,10))
  end function get_num_global_columns_f90

  function get_num_local_elems_f90 () result(num_elems) bind(c)
    use dimensions_mod,    only: nelemd
    !
    ! Local(s)
    !
    integer (kind=c_int) :: num_elems

    ! Sanity check
    call check_grids_inited(.true.)

    num_elems = nelemd
  end function get_num_local_elems_f90

  function get_num_global_elems_f90 () result(num_elems) bind(c)
    use dimensions_mod,    only: nelem
    !
    ! Local(s)
    !
    integer (kind=c_int) :: num_elems

    ! Sanity check
    call check_grids_inited(.true.)

    num_elems = nelem
  end function get_num_global_elems_f90

  function get_np_f90 () result (np_out) bind(c)
    use dimensions_mod, only: np
    !
    ! Local(s)
    !
    integer (kind=c_int) :: np_out

    np_out = np
  end function get_np_f90

  function get_nlev_f90 () result (nlev_out) bind(c)
    use dimensions_mod, only: nlev
    !
    ! Local(s)
    !
    integer (kind=c_int) :: nlev_out

    nlev_out = nlev
  end function get_nlev_f90

end module homme_grid_mod
