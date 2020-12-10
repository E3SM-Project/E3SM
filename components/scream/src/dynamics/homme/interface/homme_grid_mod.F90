module homme_grid_mod
  use iso_c_binding,     only: c_ptr, c_f_pointer, c_int, c_bool, c_double
  use parallel_mod,      only: abortmp
  use kinds,             only: iulog
  use homme_context_mod, only: par, masterproc
  use phys_grid_mod,     only: fv_nphys

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

  subroutine init_grids_f90 (pgN) bind(c)
    use dyn_grid_mod,      only: dyn_grid_init
    use phys_grid_mod,     only: phys_grid_init
    use homme_context_mod, only: is_geometry_inited
    !
    ! Input(s)
    !
    integer (kind=c_int), intent(in) :: pgN

    ! We don't expect this to be called twice
    call check_grids_inited(.false.)

    ! Init SE grid struff in Homme
    call dyn_grid_init ()

    ! If doing a dyn-only test, we don't need pg stuff
    if (pgN .ge. 0) then
      call phys_grid_init (pgN)
    endif

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

  subroutine get_dyn_grid_data_f90 (gids_ptr, elgpgp_ptr, lat_ptr, lon_ptr) bind(c)
    use homme_context_mod, only: elem
    use dimensions_mod,    only: nelemd, np
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: gids_ptr, elgpgp_ptr, lat_ptr, lon_ptr
    !
    ! Local(s)
    !
    real(kind=c_double), pointer :: lat (:), lon(:)
    integer(kind=c_int), pointer :: gids (:), elgpgp(:,:)
    integer :: idof, ip,jp, ie

    ! Sanity check
    call check_grids_inited(.true.)

    call c_f_pointer (gids_ptr,   gids,   [nelemd*np*np])
    call c_f_pointer (elgpgp_ptr, elgpgp, [3,nelemd*np*np])
    call c_f_pointer (lat_ptr,    lat,    [nelemd*np*np])
    call c_f_pointer (lon_ptr,    lon,    [nelemd*np*np])

    ! Get the gids
    idof = 1
    do ie=1,nelemd
      do jp=1,4
        do ip=1,4
          elgpgp(1,idof) = ie-1
          elgpgp(2,idof) = jp-1
          elgpgp(3,idof) = ip-1

          gids(idof) = INT(elem(ie)%gdofP(ip,jp),kind=c_int)
          lat(idof)  = elem(ie)%spherep(ip,jp)%lat
          lon(idof)  = elem(ie)%spherep(ip,jp)%lon

          idof = idof + 1
        enddo
      enddo
    enddo

  end subroutine get_dyn_grid_data_f90

  subroutine get_phys_grid_data_f90 (pg_type, gids_ptr, lat_ptr, lon_ptr, area_ptr) bind(c)
    use parallel_mod,  only: abortmp
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

    ncols = get_num_local_columns_f90()

    call c_f_pointer (gids_ptr, gids, [ncols])
    call c_f_pointer (lat_ptr,  lat,  [ncols])
    call c_f_pointer (lon_ptr,  lon,  [ncols])
    call c_f_pointer (area_ptr, area, [ncols])

    ! Retrieve all data
    call get_my_phys_data (gids, lat, lon, area, pg_type)

  end subroutine get_phys_grid_data_f90

  function get_num_local_columns_f90 () result (ncols) bind(c)
    use phys_grid_mod,     only: get_num_local_columns
    !
    ! Local(s)
    !
    integer (kind=c_int) :: ncols

    ! Sanity check
    call check_grids_inited(.true.)

    ncols = get_num_local_columns()
  end function get_num_local_columns_f90

  function get_num_global_columns_f90 () result (ncols) bind(c)
    use phys_grid_mod, only: get_num_global_columns
    !
    ! Local(s)
    !
    integer (kind=c_int) :: ncols

    ! Sanity check
    call check_grids_inited(.true.)

    ncols = get_num_global_columns()
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
