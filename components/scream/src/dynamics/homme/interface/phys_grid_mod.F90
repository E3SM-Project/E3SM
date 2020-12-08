module phys_grid_mod

  use iso_c_binding, only: c_int, c_double
  use parallel_mod,  only: abortmp, MPIinteger_t, MPIreal_t
  use kinds,         only: iulog

  implicit none
  private

  public :: phys_grid_init, cleanup_grid_init_data
  public :: finalize_phys_grid
  public :: get_my_phys_data
  public :: get_num_local_columns, get_num_global_columns

  ! Available options for pg balancing
  integer, public, parameter :: pg_default   = 0
  integer, public, parameter :: pg_twin_cols = 1  ! Not yet supported in SCREAM!

  ! Global geo info
  real (kind=c_double), pointer :: g_lat(:)
  real (kind=c_double), pointer :: g_lon(:)
  real (kind=c_double), pointer :: g_area(:)
  integer (kind=c_int), pointer :: g_dofs(:)
  integer (kind=c_int), pointer :: g_ncols(:)
  integer (kind=c_int), pointer :: g_offsets(:)

  integer, public :: fv_nphys = 0

  logical :: is_phys_grid_inited = .false.

contains

  subroutine check_phys_grid_inited ()
    if (.not. is_phys_grid_inited) then
      call abortmp ("Error! Physics grid was not inited.\n")
    endif
  end subroutine check_phys_grid_inited

  subroutine phys_grid_init (pgN)
    use gllfvremap_mod,    only: gfr_init
    use homme_context_mod, only: elem, par
    !
    ! Input(s)
    !
    integer (kind=c_int), intent(in) :: pgN
    !
    ! Local(s)
    !
    integer :: ldofs, ierr, proc, ngcols

    if (pgN>0) then
      fv_nphys = pgN
      call gfr_init(par, elem, fv_nphys)
    endif

    ! Set this right away, so calls to get_num_[local|global]_columns doesn't crap out
    is_phys_grid_inited = .true.

    ngcols = get_num_global_columns ()

    allocate(g_area(ngcols))
    allocate(g_lat(ngcols))
    allocate(g_lon(ngcols))
    allocate(g_dofs(ngcols))
    allocate(g_ncols(par%nprocs))
    allocate(g_offsets(par%nprocs))

    ! Gather local counts, then do a scan sum
    ldofs = get_num_local_columns ()
    g_ncols(:) = 0
    call MPI_Allgather(ldofs, 1, MPIinteger_t, g_ncols, 1, MPIinteger_t, par%comm, ierr)
    g_offsets(:) = 0
    do proc=2,par%nprocs,1
      g_offsets(proc) = g_offsets(proc-1) + g_ncols(proc-1)
    enddo

    ! Fill arrays of global coords, area, dofs
    call compute_global_coords (g_lat, g_lon)
    call compute_global_area (g_area)
    call compute_global_dofs (g_dofs)

  end subroutine phys_grid_init

  subroutine cleanup_grid_init_data ()
    if (is_phys_grid_inited) then
      deallocate(g_ncols)
    endif
  end subroutine cleanup_grid_init_data

  subroutine finalize_phys_grid ()
    use gllfvremap_mod, only: gfr_finish

    if (fv_nphys>0) then
      call gfr_finish()
    endif

    deallocate(g_area)
    deallocate(g_lat)
    deallocate(g_lon)
    deallocate(g_dofs)
    deallocate(g_offsets)
  
    is_phys_grid_inited = .false.
  end subroutine finalize_phys_grid

  function get_num_local_columns () result (ncols)
    use dimensions_mod,    only: nelemd
    use homme_context_mod, only: elem
    !
    ! Local(s)
    !
    integer :: ncols, ie

    ! Sanity check
    call check_phys_grid_inited()

    if (fv_nphys>0) then
      ncols = nelemd*fv_nphys*fv_nphys
    else
      ncols = 0
      do ie=1,nelemd
        ncols = ncols + elem(ie)%idxP%NumUniquePts
      enddo
    endif
  end function get_num_local_columns

  function get_num_global_columns () result (num_cols) bind(c)
    use homme_context_mod, only: is_geometry_inited, elem
    use dimensions_mod,    only: nelem, np
    !
    ! Local(s)
    !
    integer (kind=c_int) :: num_cols
    integer :: ie

    ! Sanity check
    call check_phys_grid_inited()

    if (fv_nphys>0) then
      num_cols = nelem*fv_nphys*fv_nphys
    else
      num_cols = nelem*(np-1)*(np-1)+2
    endif
  end function get_num_global_columns

  subroutine get_my_phys_data (gids, lat, lon, area, pg_type)
    use homme_context_mod, only: elem, iam
    use dimensions_mod,    only: nelemd
    !
    ! Input(s)
    !
    integer(kind=c_int), pointer, intent(in) :: gids(:)
    real(kind=c_double), pointer, intent(in) :: lat(:), lon(:), area(:)
    integer(kind=c_int),          intent(in) :: pg_type
    !
    ! Local(s)
    !
    integer :: pgN, balancing_opt, idof, ndofs

    call check_phys_grid_inited ()

    ! Possible values for pg_type:
    !   0 : physics grid on GLL nodes
    !   2 : physics grid on FV points in 2x2 subcell (PG2)
    !   3 : physics grid on FV points in 3x3 subcell (PG3)
    !   4 : physics grid on FV points in 4x4 subcell (PG4)
    !  10 : GLL points, use twin columns
    !  12 : PG2 points, use twin columns
    !  13 : PG3 points, use twin columns
    !  14 : PG4 points, use twin columns
    ! In general:
    !   *0 => GLL nodes, *N,N>0 => FV points
    !   1* => redistribute using twin columns

    pgN = mod(pg_type,10)
    if (pgN .ne. fv_nphys) then
      call abortmp ("Error! Requested N for phys grid different from what was used at init time.")
    endif

    balancing_opt = pg_type / 10

    if (balancing_opt .ne. pg_default) then
      call abortmp ("Error! Twin columns not yet implemented for physics grid.")
    else
      ! TODO: when you enable twin columns, you'll have to manually
      !       do the search, since you can't just grab the offset-ed entries
      ndofs = get_num_local_columns ()
      do idof=1,ndofs
        gids(idof) = g_dofs(g_offsets(iam+1) + idof)
        lat(idof)  = g_lat(g_offsets(iam+1) + idof)
        lon(idof)  = g_lon(g_offsets(iam+1) + idof)
        area(idof) = g_area(g_offsets(iam+1) + idof)
      enddo

    endif
  end subroutine get_my_phys_data

  subroutine compute_global_dofs (dofs_g)
    use dimensions_mod,    only: nelemd
    use homme_context_mod, only: elem, par, iam
    !
    ! Input(s)
    !
    integer (kind=c_int), pointer, intent(inout) :: dofs_g(:)
    !
    ! Local(s)
    !
    integer (kind=c_int), pointer :: dofs_l(:)
    integer :: ie, offset, i, j, icol, idof, ncols, ierr

    if (fv_nphys .gt. 0) then
      ! TODO
      call abortmp ("Error! Phys grid N>0 not yet implemented in scream.")
    else
      offset = g_offsets(iam+1)

      ncols  = g_ncols(iam+1)
      dofs_l => dofs_g(offset+1:offset+ncols)

      ! Fill local dofs part
      idof = 1
      do ie=1,nelemd
        do icol=1,elem(ie)%idxP%NumUniquePts
          i = elem(ie)%idxP%ia(icol)
          j = elem(ie)%idxP%ja(icol)

          dofs_l(idof) = elem(ie)%gdofP(i,j)
          idof = idof + 1
        enddo
      enddo

      ! Gather all the other dofs
      call MPI_Allgatherv( dofs_l, ncols, MPIinteger_t, &
                           dofs_g, g_ncols, g_offsets, MPIinteger_t, par%comm, ierr)
    endif

  end subroutine compute_global_dofs

  subroutine compute_global_area(area_g)
    use dof_mod,           only: UniquePoints
    use dimensions_mod,    only: np, nelemd
    use homme_context_mod, only: elem, par, iam, masterproc
    !
    ! Input(s)
    !
    real(kind=c_double), pointer, intent(in) :: area_g(:)
    !
    ! Local(s)
    !
    real(kind=c_double), pointer :: area_l(:)
    real(kind=c_double), dimension(np,np)  :: areaw
    integer  :: ie, offset, start, ierr, ncols

    if (masterproc) then
      write(iulog,*) 'INFO: Non-scalable action: Computing global area in SE dycore.'
    endif

    if (fv_nphys > 0) then
      ! TODO
      call abortmp("Error! Not yet implemented. Why didn't you get an error before?")
    else
      ! physics is on GLL grid
      offset = g_offsets(iam+1)

      ncols = get_num_local_columns()
      area_l => area_g(offset+1 : offset+ncols)
      do ie=1,nelemd
        areaw = 1.0_c_double / elem(ie)%rspheremp(:,:)
        ncols = elem(ie)%idxP%NumUniquePts
        start = elem(ie)%idxP%UniquePtOffset - offset
        call UniquePoints(elem(ie)%idxP, areaw, area_l(start:start+ncols-1))
      enddo

      ncols = get_num_local_columns()
      call MPI_Allgatherv(area_l, ncols, MPIreal_t, &
                          area_g, g_ncols, g_offsets, MPIreal_t, &
                          par%comm, ierr)
    endif

  end subroutine compute_global_area

  subroutine compute_global_coords(lat_g, lon_g)
    use dimensions_mod,    only: nelemd
    use dof_mod,           only: UniqueCoords
    use gllfvremap_mod,    only: gfr_f_get_latlon
    use homme_context_mod, only: elem, par, iam, masterproc
    !
    ! Input(s)
    !
    real(kind=c_double), pointer, intent(out) :: lat_g(:)
    real(kind=c_double), pointer, intent(out) :: lon_g(:)
    !
    ! Local(s)
    !
    real(kind=c_double), pointer :: lat_l(:), lon_l(:)
    integer  :: ncols, ie, offset, start, ierr

    if (masterproc) then
      write(iulog,*) 'INFO: Non-scalable action: Computing global coords in SE dycore.'
    end if

    if (fv_nphys > 0) then
      ! TODO
      call abortmp("Error! Not yet implemented. Why didn't you get an error before?")
    else
      ! physics is on GLL grid

      offset = g_offsets(iam+1)

      ncols = get_num_local_columns()
      lat_l => lat_g(offset+1 : offset+ncols)
      lon_l => lon_g(offset+1 : offset+ncols)

      do ie=1,nelemd
        ncols = elem(ie)%idxP%NumUniquePts
        start = elem(ie)%idxP%UniquePtOffset - offset
        call UniqueCoords(elem(ie)%idxP, elem(ie)%spherep, &
                          lat_l(start:start+ncols-1),      &
                          lon_l(start:start+ncols-1))
      enddo

      ncols = get_num_local_columns()

      call MPI_Allgatherv(lat_l, ncols,   MPIreal_t, &
                          lat_g, g_ncols, g_offsets, MPIreal_t, &
                          par%comm, ierr)
      call MPI_Allgatherv(lon_l, ncols, MPIreal_t, &
                          lon_g, g_ncols, g_offsets, MPIreal_t, &
                          par%comm, ierr)
    end if

  end subroutine compute_global_coords

end module phys_grid_mod
