module phys_grid_mod

  use iso_c_binding, only: c_int, c_double
  use parallel_mod,  only: abortmp, MPIinteger_t
  use kinds,         only: iulog
  use mpi

  implicit none
  private

  public :: phys_grids_init, cleanup_grid_init_data
  public :: finalize_phys_grid
  public :: get_my_phys_data
  public :: get_num_local_columns, get_num_global_columns

  ! Available options for pg balancing
  integer, public, parameter :: load_bal_none      = 0
  integer, public, parameter :: load_bal_twin_cols = 1  ! Not yet supported in SCREAM!

  ! Max phys grid N. Convention: pgN=0 -> GLL grid, pgN=1,2,3,4 -> FV grid
  integer, public, parameter :: pgN_min = 0
  integer, public, parameter :: pgN_max = 4

  ! Global grid info
  ! Arrays of size nprocs
  integer (kind=c_int), pointer :: g_elem_per_rank(:)     ! Number of elems on each rank
  integer (kind=c_int), pointer :: g_elem_offsets(:)      ! Offset of each rank in global nelem-sized arrays
  ! Arrays of size nelem
  integer (kind=c_int), pointer :: g_elem_gids(:)         ! List of elem gids, spliced by rank

  ! This struct holds the info of a physics grid.
  ! NOTE: there is one struct for each pg type,
  type :: pg_specs_t
    ! Arrays of size nprocs
    integer (kind=c_int), pointer :: g_dofs_per_rank(:)     ! Number of cols on each rank
    integer (kind=c_int), pointer :: g_dofs_offsets(:)      ! Per-rank column offset; retained for potential future rebalancing support

    ! Local arrays (size = number of columns on this rank only).
    ! Each rank stores only its own columns, avoiding the large global allocation.
    real (kind=c_double), pointer :: l_lat(:)     ! Latitude coordinates
    real (kind=c_double), pointer :: l_lon(:)     ! Longitude coordinates
    real (kind=c_double), pointer :: l_area(:)    ! Column area
    integer (kind=c_int), pointer :: l_dofs(:)    ! Column global index

    ! Arrays of size nelem
    integer (kind=c_int), pointer :: g_dofs_per_elem(:)     ! List of num dofs on each elem, spliced by rank

    integer :: pgN
    logical :: inited = .false.
  end type pg_specs_t

  type(pg_specs_t), target :: pg_specs (pgN_min:pgN_max)

contains

  subroutine check_phys_grid_inited (pgN)
    !
    ! Inputs(s)
    !
    integer, intent(in) :: pgN  ! The phys grid type to check
    !
    ! Local(s)
    !
    character(1) :: str

    write(str, '(I1)') pgN
    if (pgN .lt. pgN_min .or. pgN .gt. pgN_max) then
      call abortmp ("Error! Physics grid '"//str//"' is not supported.\n")
    endif

    if (.not. pg_specs(pgN)%inited) then
      call abortmp ("Error! Physics grid '"//str//"' was not inited.\n")
    endif
  end subroutine check_phys_grid_inited

  subroutine phys_grids_init (pg_types)
    use gllfvremap_mod,    only: gfr_init
    use homme_context_mod, only: elem, par
    use dimensions_mod,    only: nelem, nelemd
    !
    ! Input(s)
    !
    integer (kind=c_int), intent(in) :: pg_types(:)
    !
    ! Local(s)
    !
    integer :: ie, ierr, proc, i, load_bal, pgN, fvN
    logical :: do_gll
    character(2) :: str

    fvN = -1
    do_gll = .false.
    do i=1, size(pg_types)
      pgN = mod(pg_types(i),10)
      load_bal = pg_types(i) / 10

      if (load_bal .ne. load_bal_none) then
        write(str,'(I2)') load_bal
        call abortmp("Error! Load balancing option '"//str//"' not supported.")
      endif

      if (pgN .lt. pgN_min .or. pgN .gt. pgN_max) then
        write(str,'(I2)') pgN
        call abortmp("Error! Physics grid type '"//str//"' not supported.")
      endif

      if (pgN > 0) then
        if (fvN .ne. -1 .and. fvN .ne. pgN) then
          call abortmp("Error! Only ONE finite volume phys grid can be enabled.")
        endif
        fvN = pgN
      else
        do_gll = .true.
      endif
    enddo
    
    ! Compute elem-related quantities, which are common to all phys grids

    ! Gather num elems on each rank
    allocate(g_elem_per_rank(par%nprocs))
    call MPI_Allgather(nelemd, 1, MPIinteger_t, g_elem_per_rank, 1, MPIinteger_t, par%comm, ierr)

    ! Compute offset of each rank in the g_elem_per_rank array
    ! This allows each rank to know where to fill arrays of size num_global_elements
    allocate(g_elem_offsets(par%nprocs))
    g_elem_offsets = 0
    do proc=2,par%nprocs
      g_elem_offsets(proc) = g_elem_offsets(proc-1) + g_elem_per_rank(proc-1)
    enddo

    ! Gather elem GIDs on each rank
    allocate(g_elem_gids(nelem))
    do ie=1,nelemd
      g_elem_gids(g_elem_offsets(par%rank+1)+ie) = elem(ie)%globalID
    enddo
    call MPI_Allgatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                         g_elem_gids, g_elem_per_rank, g_elem_offsets, MPIinteger_t, par%comm, ierr)

    ! Now init grid-specific quantities for the requested pg types
    if (do_gll) then
      call phys_grid_init (0)
    endif
    if (fvN .gt. 0) then
      call phys_grid_init (fvN)
    endif
  end subroutine phys_grids_init

  subroutine cleanup_grid_init_data ()
    !
    ! Local(s)
    !
    integer :: pgN

    ! Cleanup data no longer needed in each pg_specs_t struct
    do pgN = pgN_min,pgN_max
      if (pg_specs(pgN)%inited) then
        ! nprocs-sized arrays
        deallocate(pg_specs(pgN)%g_dofs_per_rank)

        ! nelem-sized arrays
        deallocate(pg_specs(pgN)%g_dofs_per_elem)
      endif
    enddo

    ! Cleanup global elem-related arrays
    deallocate(g_elem_per_rank)
    deallocate(g_elem_offsets)
    deallocate(g_elem_gids)
  end subroutine cleanup_grid_init_data

  subroutine finalize_phys_grid ()
    use gllfvremap_mod, only: gfr_finish
    !
    ! Local(s)
    !
    integer :: pgN

    do pgN = pgN_min,pgN_max
      if (pg_specs(pgN)%inited) then
        deallocate(pg_specs(pgN)%l_area)
        deallocate(pg_specs(pgN)%l_lat)
        deallocate(pg_specs(pgN)%l_lon)
        deallocate(pg_specs(pgN)%l_dofs)

        deallocate(pg_specs(pgN)%g_dofs_offsets)

        if (pgN>0) then
          call gfr_finish()
        endif

        pg_specs(pgN)%inited = .false.
      endif
    enddo
  end subroutine finalize_phys_grid

  function get_num_local_columns (pgN) result (ncols)
    use dimensions_mod,    only: nelemd
    use homme_context_mod, only: elem
    !
    ! Input(s)
    !
    integer, intent(in) :: pgN
    !
    ! Local(s)
    !
    integer :: ncols, ie

    ! Sanity check
    call check_phys_grid_inited(pgN)

    if (pgN>0) then
      ncols = nelemd*pgN*pgN
    else
      ncols = 0
      do ie=1,nelemd
        ncols = ncols + elem(ie)%idxP%NumUniquePts
      enddo
    endif
  end function get_num_local_columns

  function get_num_global_columns (pgN) result (num_cols) bind(c)
    use dimensions_mod,    only: nelem, np
    !
    ! Input(s)
    !
    integer (kind=c_int), intent(in) :: pgN
    !
    ! Local(s)
    !
    integer (kind=c_int) :: num_cols

    ! Sanity check
    call check_phys_grid_inited(pgN)

    if (pgN>0) then
      num_cols = nelem*pgN*pgN
    else
      num_cols = nelem*(np-1)*(np-1)+2
    endif
  end function get_num_global_columns

  subroutine get_my_phys_data (gids, lat, lon, area, pg_type)
    use shr_const_mod,     only: pi=>SHR_CONST_PI
    use control_mod,       only: geometry
    !
    ! Input(s)
    !
    integer(kind=c_int), intent(out) :: gids(:)
    real(kind=c_double), intent(out) :: lat(:), lon(:), area(:)
    integer(kind=c_int), intent(in) :: pg_type
    !
    ! Local(s)
    !
    integer :: pgN, load_bal, idof, ndofs
    type(pg_specs_t), pointer :: pgs
    logical :: is_sphere

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

    call check_phys_grid_inited (pgN)

    pgs => pg_specs(pgN)

    load_bal = pg_type / 10

    if (load_bal .ne. load_bal_none) then
      call abortmp ("Error! No load balancing option implemented in scream yet.")
    else
      ndofs = get_num_local_columns (pgN)
      is_sphere = trim(geometry) /= 'plane'
      do idof=1,ndofs
        gids(idof)  = pgs%l_dofs(idof)
        lat(idof)   = pgs%l_lat(idof)
        lon(idof)   = pgs%l_lon(idof)
        if (is_sphere) then
          lat(idof) = lat(idof) * 180.0_c_double / pi
          lon(idof) = lon(idof) * 180.0_c_double / pi
        end if
        area(idof)  = pgs%l_area(idof)
      enddo

    endif
  end subroutine get_my_phys_data

  subroutine compute_local_dofs (pg)
    use dimensions_mod,    only: nelem
    use homme_context_mod, only: par
    !
    ! Input(s)
    !
    type(pg_specs_t), intent(inout) :: pg
    !
    ! Local(s)
    !
    integer :: ie, icol, idof, proc, elem_gid, elem_offset
    integer, allocatable :: exclusive_scan_dofs_per_elem(:)

    ! Compute the prefix sum of dofs per element (indexed by global elem ID).
    ! This gives globally-unique column IDs even though we only store local data.
    allocate(exclusive_scan_dofs_per_elem(nelem))
    exclusive_scan_dofs_per_elem(1) = 0
    do ie=2,nelem
      exclusive_scan_dofs_per_elem(ie) = exclusive_scan_dofs_per_elem(ie-1) + pg%g_dofs_per_elem(ie-1)
    enddo

    ! Build the per-rank column offset array. Not used for local-only data access,
    ! but retained for potential future rebalancing support (e.g., twin columns).
    allocate (pg%g_dofs_offsets(par%nprocs))
    pg%g_dofs_offsets(1)=0
    do proc=2,par%nprocs
      pg%g_dofs_offsets(proc) = pg%g_dofs_offsets(proc-1) + pg%g_dofs_per_rank(proc-1)
    enddo

    ! Allocate local-only dofs array and fill with globally-unique column IDs
    ! for elements owned by this rank.
    allocate(pg%l_dofs(get_num_local_columns(pg%pgN)))
    idof = 1
    elem_offset = g_elem_offsets(par%rank+1)
    do ie=1,g_elem_per_rank(par%rank+1)
      elem_gid = g_elem_gids(elem_offset+ie)
      do icol=1,pg%g_dofs_per_elem(elem_gid)
        pg%l_dofs(idof) = exclusive_scan_dofs_per_elem(elem_gid)+icol
        idof = idof+1
      enddo
    enddo
  end subroutine compute_local_dofs

  subroutine compute_local_area(pg)
    use dof_mod,           only: UniquePoints
    use dimensions_mod,    only: np, nelemd
    use homme_context_mod, only: elem
    use gllfvremap_mod,    only: gfr_f_get_area
    !
    ! Input(s)
    !
    type(pg_specs_t), intent(inout) :: pg
    !
    ! Local(s)
    !
    real(kind=c_double), dimension(np,np)  :: areaw
    integer :: ie, start, ncols, i, j, k

    allocate(pg%l_area(get_num_local_columns(pg%pgN)))

    if (pg%pgN > 0) then
       start = 0
       do ie=1,nelemd
          do j = 1,pg%pgN
             do i = 1,pg%pgN
                k = start + (j-1)*pg%pgN + i
                pg%l_area(k) = gfr_f_get_area(ie, i, j)
             end do
          end do
          start = start + pg%pgN*pg%pgN
       end do
    else
      ! physics is on GLL grid
      start = 1
      do ie=1,nelemd
        areaw = 1.0_c_double / elem(ie)%rspheremp(:,:)
        ncols = elem(ie)%idxP%NumUniquePts
        call UniquePoints(elem(ie)%idxP, areaw, pg%l_area(start:start+ncols-1))
        start = start + ncols
      enddo
    endif
  end subroutine compute_local_area

  subroutine compute_local_coords(pg)
    use dimensions_mod,    only: nelemd
    use dof_mod,           only: UniqueCoords
    use homme_context_mod, only: elem
    use gllfvremap_mod,    only: gfr_f_get_latlon
    !
    ! Input(s)
    !
    type(pg_specs_t), intent(inout) :: pg
    !
    ! Local(s)
    !
    integer  :: ncols, ie, start, i, j, k

    allocate(pg%l_lat(get_num_local_columns(pg%pgN)))
    allocate(pg%l_lon(get_num_local_columns(pg%pgN)))

    if (pg%pgN > 0) then
       start = 0
       do ie=1,nelemd
          do j = 1,pg%pgN
             do i = 1,pg%pgN
                k = start + (j-1)*pg%pgN + i
                call gfr_f_get_latlon(ie, i, j, pg%l_lat(k), pg%l_lon(k))
             end do
          end do
          start = start + pg%pgN*pg%pgN
       end do
    else
      ! physics is on GLL grid
      start = 1
      do ie=1,nelemd
        ncols = elem(ie)%idxP%NumUniquePts
        call UniqueCoords(elem(ie)%idxP, elem(ie)%spherep, &
                          pg%l_lat(start:start+ncols-1),   &
                          pg%l_lon(start:start+ncols-1))
        start = start + ncols
      enddo
    endif
  end subroutine compute_local_coords

  subroutine phys_grid_init (pgN)
    use gllfvremap_mod,    only: gfr_init
    use homme_context_mod, only: elem, par
    use dimensions_mod,    only: nelem, nelemd
#ifdef HAVE_MOAB
    use seq_comm_mct,      only: MHID, MHFID  ! id of homme moab coarse and fine applications
    use seq_comm_mct,      only: ATMID
    use seq_comm_mct,      only: mhpgid       ! id of pgx moab application
    use semoab_mod,        only: create_moab_meshes
    use iMOAB, only : iMOAB_RegisterApplication
    use iso_c_binding
#endif
    !
    ! Input(s)
    !
    integer, intent(in), target :: pgN
    !
    ! Local(s)
    !
    integer, allocatable :: dofs_per_elem (:)
    integer :: ierr, proc, ie, offset
    character(2) :: str
    type(pg_specs_t), pointer :: pg

#ifdef HAVE_MOAB
    integer :: ATM_ID1
    character*32  appname
#endif
    pg => pg_specs(pgN)

    if (pg%inited) then
      write (str, '(I1)') pgN
      call abortmp("Error! Physics grid with pgN="//str//" was already inited.")
    endif

    ! Set this to true immediately, so we can call get_num_xyz_columns
    pg%inited = .true.
    pg%pgN = pgN

    if (pgN>0) call gfr_init(par, elem, pgN)

    ! Gather the number of unique cols on each rank
    allocate(pg%g_dofs_per_rank(par%nprocs))
    call MPI_Allgather(get_num_local_columns(pg%pgN), 1, MPIinteger_t, &
                       pg%g_dofs_per_rank, 1, MPIinteger_t, par%comm, ierr)

    ! Gather the number of unique cols on each element
    ! NOTE: we use a temp (dofs_per_elem) rather than g_dofs, since in g_dofs we order by
    !       elem GID. But elem GIDs may not be distributed linearly across ranks
    !       (e.g., proc0 may own [1-20,41-60],and proc1 may owns [21-40]), while in order
    !       to use Allgatherv, we need the send buffer to be contiguous. So in dofs_per_elem
    !       we order the data *by rank*. Once the gather is complete, we copy the data from
    !       dofs_per_elem into g_dofs_per_elem, ordering by elem GID instead.
    !
    allocate(dofs_per_elem(nelem))
    if (pgN>0) then
       do ie=1,nelem
          dofs_per_elem(ie) = pgN*pgN
       end do
    else
       do ie=1,nelemd
          dofs_per_elem(g_elem_offsets(par%rank+1)+ie) = elem(ie)%idxP%NumUniquePts
       enddo
       call MPI_Allgatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
            dofs_per_elem, g_elem_per_rank, g_elem_offsets, MPIinteger_t, par%comm, ierr)
    end if

    allocate(pg%g_dofs_per_elem(nelem))
    do proc=1,par%nprocs
      offset = g_elem_offsets(proc)
      do ie=1,g_elem_per_rank(proc)
        pg%g_dofs_per_elem(g_elem_gids(offset+ie)) = dofs_per_elem(offset+ie)
      enddo
    enddo

    ! Compute local arrays of coords, area, and globally-unique dof IDs
    call compute_local_dofs (pg)
    call compute_local_coords (pg)
    call compute_local_area (pg)
#ifdef HAVE_MOAB
    if (pgN > 0) then
        appname="HM_COARSE"//C_NULL_CHAR
        ATM_ID1 = 120 !
        ierr = iMOAB_RegisterApplication(appname, par%comm, ATM_ID1, MHID)
        if (ierr > 0 )  &
        call abortmp('Error: cannot register moab app')
        if(par%masterproc) then
            write(iulog,*) " "
            write(iulog,*) "register MOAB app:", trim(appname), "  MHID=", MHID
            write(iulog,*) " "
        endif
        appname="HM_FINE"//C_NULL_CHAR
        ATM_ID1 = 119 ! this number should not conflict with other components IDs; how do we know?
        ierr = iMOAB_RegisterApplication(appname, par%comm, ATM_ID1, MHFID)
        if (ierr > 0 )  &
        call abortmp('Error: cannot register moab app for fine mesh')
        if(par%masterproc) then
            write(iulog,*) " "
            write(iulog,*) "register MOAB app:", trim(appname), "  MHFID=", MHFID
            write(iulog,*) " "
        endif
        appname="HM_PGX"//C_NULL_CHAR
        ATM_ID1 =  ATMID(1) ! this number should not conflict with other components IDs; how do we know?
        !
        ! in this case, we reuse the main atm id, mhid will not be used for intersection anymore
        ! still, need to be careful
        ierr = iMOAB_RegisterApplication(appname, par%comm, ATM_ID1, mhpgid)
        if (ierr > 0 )  &
            call abortmp('Error: cannot register moab app for fine mesh')
        if(par%masterproc) then
            write(iulog,*) " "
            write(iulog,*) "register MOAB app:", trim(appname), "  MHPGID=", mhpgid
            write(iulog,*) " "
        endif
! instance distributed moab meshes from elem structures
!    1 ) spectral coarse mesh
!    2 ) GLL fine quad mesh (used mostly for visualization)
!    3 ) pgN FV type mesh, (most of the time pg2 mesh), used for coupling with other components;
        call create_moab_meshes(par, elem, pgN)
    endif
#endif
  end subroutine phys_grid_init

end module phys_grid_mod
