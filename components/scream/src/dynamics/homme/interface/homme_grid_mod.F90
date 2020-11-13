module homme_grid_mod
  use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_bool
  use parallel_mod,  only: abortmp

  implicit none
  private

  public :: init_geometry_f90

  public :: get_cols_gids_f90, get_cols_specs_f90
  public :: get_num_owned_columns_f90
  public :: get_num_owned_elems_f90, get_np_f90, get_nlev_f90

contains

  subroutine init_geometry_f90 () bind(c)
    use prim_driver_base,     only: prim_init1_geometry, prim_init1_cleanup, MetaVertex, GridEdge
    use prim_cxx_driver_base, only: init_cxx_connectivity
    use dimensions_mod,       only: nelemd
    use homme_context_mod,    only: is_geometry_inited, is_parallel_inited, &
                                    elem, par, dom_mt

    if (.not. is_parallel_inited) then
      call abortmp ("Error! 'homme_init_parallel_f90' must be called *before* init_geometry_f90.\n")
    endif

    print *, "Initing geometry..."

    call prim_init1_geometry(elem,par,dom_mt)

    call init_cxx_connectivity(nelemd,GridEdge,MetaVertex,par)

    ! Cleanup the tmp stuff used in prim_init1_geometry
    call prim_init1_cleanup()

    is_geometry_inited = .true.
  end subroutine init_geometry_f90

  subroutine get_cols_gids_f90 (gids_ptr, owned_only) bind(c)
    use homme_context_mod, only: is_geometry_inited
    use dimensions_mod,    only: np, nelemd
    use element_mod,       only: index_t
    use kinds,             only: long_kind
    use homme_context_mod, only: elem

    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: gids_ptr
    logical (kind=c_bool), intent(in) :: owned_only
    !
    ! Local(s)
    !
    type (index_t) :: idxP
    integer (kind=long_kind), pointer :: gids(:)
    integer :: i,ip,jp,ie,ncols,icol

    if (.not. is_geometry_inited) then
      call abortmp ("Error! Geometry was not inited yet.\n")
    endif

    if (owned_only) then
      ncols = get_num_owned_columns_f90()
    else
      ncols = nelemd*np*np
    endif

    call c_f_pointer(gids_ptr, gids, [ncols])

    icol = 1
    do ie=1,nelemd
      idxP = elem(ie)%idxP
      if (owned_only) then
        ncols = idxP%NumUniquePts
      else
        ncols = np*np
      endif
      do i=1,ncols
        ip = idxP%ia(i)
        jp = idxP%ja(i)

        gids(icol) = elem(ie)%gdofP(ip,jp)-1
        icol = icol+1
      enddo
    enddo
  end subroutine get_cols_gids_f90

  subroutine get_cols_specs_f90 (gids_ptr, elgp_ptr, owned_only) bind(c)
    use homme_context_mod, only: is_geometry_inited, elem
    use dimensions_mod,    only: np, nelemd
    use element_mod,       only: index_t
    use kinds,             only: long_kind
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: gids_ptr, elgp_ptr
    logical (kind=c_bool), intent(in) :: owned_only
    !
    ! Local(s)
    !
    type (index_t) :: idxP
    integer (kind=long_kind), pointer :: gids(:)
    integer (kind=c_int), pointer :: elgp(:,:)
    integer :: i,ip,jp,ie,ncols,icol

    if (.not. is_geometry_inited) then
      call abortmp ("Error! Geometry was not inited yet.\n")
    endif

    if (owned_only) then
      ncols = get_num_owned_columns_f90()
    else
      ncols = nelemd*np*np
    endif

    call c_f_pointer(gids_ptr, gids, [ncols])
    call c_f_pointer(elgp_ptr, elgp, [3,ncols])

    icol = 1
    do ie=1,nelemd
      idxP = elem(ie)%idxP
      if (owned_only) then
        ncols = idxP%NumUniquePts
        do i=1,ncols
          ip = idxP%ia(i)
          jp = idxP%ja(i)
          elgp(1,icol) = ie-1
          elgp(2,icol) = jp-1
          elgp(3,icol) = ip-1

          gids(icol) = elem(ie)%gdofP(ip,jp)-1
          icol = icol + 1
        enddo
      else
        do jp=1,4
          do ip=1,4
            elgp(1,icol) = ie-1
            elgp(2,icol) = jp-1
            elgp(3,icol) = ip-1

            gids(icol) = elem(ie)%gdofP(ip,jp)-1
            icol = icol + 1
          enddo
        enddo
      endif
    enddo
  end subroutine get_cols_specs_f90

  function get_num_owned_columns_f90 () result (num_cols) bind(c)
    use homme_context_mod, only: is_geometry_inited, elem
    use dimensions_mod,    only: nelemd
    !
    ! Local(s)
    !
    integer (kind=c_int) :: num_cols
    integer :: ie

    if (.not. is_geometry_inited) then
      call abortmp ("Error! Geometry was not inited yet.\n")
    endif

    num_cols = 0
    do ie=1,nelemd
      num_cols = num_cols + elem(ie)%idxP%NumUniquePts
    enddo
  end function get_num_owned_columns_f90

  function get_num_owned_elems_f90 () result(num_elems) bind(c)
    use homme_context_mod, only: is_geometry_inited
    use dimensions_mod,    only: nelemd
    !
    ! Local(s)
    !
    integer (kind=c_int) :: num_elems

    if (.not. is_geometry_inited) then
      call abortmp ("Error! Geometry was not inited yet.\n")
    endif

    print *, "nelemd: ", nelemd
    num_elems = nelemd
  end function get_num_owned_elems_f90

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
