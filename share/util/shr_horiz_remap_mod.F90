module shr_horiz_remap_mod
  !-------------------------------------------------------------------------------------------
  !
  ! Shared horizontal remapping infrastructure for E3SM components.
  ! Inspired by the EAMxx horiz remap code...
  !
  ! Reads ESMF/TempestRemap mapping files and applies sparse matrix-vector
  ! multiply to remap fields from a source grid to a target lat-lon grid.
  ! This module contains the component-independent core algorithm; each
  ! component provides a thin wrapper that handles its own decomposition.
  !
  ! The init is split into three phases to separate shared from component
  ! logic:
  !   Phase 1: shr_horiz_remap_read_mapfile  - read map file, partition target grid
  !   Phase 2: shr_horiz_remap_build_comm    - build MPI comm pattern from ownership array
  !   Phase 3: shr_horiz_remap_apply         - runtime remap (send_buf already packed)
  !
  ! Each component wrapper:
  !   1. Calls phase 1 (shared: reads map file)
  !   2. Fills gcol_to_rank(1:n_a) using its own decomposition info
  !   3. Calls phase 2 (shared: builds MPI pattern, returns send_gcol_list)
  !   4. Converts send_gcol_list to local indices using its own indexing
  !   5. At runtime: packs send_buf from its own data layout, calls phase 3
  !
  !-------------------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private
  save

  public :: shr_horiz_remap_t
  public :: shr_horiz_remap_read_mapfile
  public :: shr_horiz_remap_build_comm
  public :: shr_horiz_remap_apply

  type shr_horiz_remap_t
    logical  :: initialized = .false.
    integer  :: n_a = 0          ! source grid size
    integer  :: n_b = 0          ! target grid size
    integer  :: nlat = 0         ! target latitude dimension
    integer  :: nlon = 0         ! target longitude dimension
    integer  :: n_b_local = 0    ! number of target points owned by this rank
    integer  :: row_start = 0    ! global starting row (1-indexed) for this rank

    ! Local sparse matrix (triplet form, only rows owned by this rank)
    integer  :: nnz_local = 0
    integer,  allocatable :: dst_local(:)  ! local target index (1..n_b_local)
    integer,  allocatable :: src_gid(:)    ! global source column ID (1-indexed)
    real(r8), allocatable :: wgt(:)        ! interpolation weight

    ! Source data gathering: src_need_gids is sorted for binary search
    integer  :: n_src_need = 0
    integer,  allocatable :: src_need_gids(:)    ! sorted unique global IDs of needed src columns
    integer,  allocatable :: src_need_recvidx(:) ! recv buffer index for each src_need_gids entry

    ! MPI communication pattern for Alltoallv
    integer,  allocatable :: send_counts(:)  ! (0:npes-1)
    integer,  allocatable :: send_displs(:)  ! (0:npes-1)
    integer,  allocatable :: recv_counts(:)  ! (0:npes-1)
    integer,  allocatable :: recv_displs(:)  ! (0:npes-1)
    integer  :: n_send_total = 0
    integer  :: n_recv_total = 0

    ! Target grid coordinates
    real(r8), allocatable :: lat(:)   ! nlat
    real(r8), allocatable :: lon(:)   ! nlon
  end type shr_horiz_remap_t

CONTAINS

  !-------------------------------------------------------------------------------------------
  pure integer function bsearch(arr, n, val)
    ! Binary search for val in sorted array arr(1:n). Returns index or 0.
    integer, intent(in) :: n, val
    integer, intent(in) :: arr(n)
    integer :: lo, hi, mid
    lo = 1; hi = n
    bsearch = 0
    do while (lo <= hi)
      mid = (lo + hi) / 2
      if (arr(mid) == val) then
        bsearch = mid
        return
      else if (arr(mid) < val) then
        lo = mid + 1
      else
        hi = mid - 1
      end if
    end do
  end function bsearch

  !-------------------------------------------------------------------------------------------
  pure integer function src_gid_to_recvidx(rd, gid)
    ! Map a global source column ID to its position in the recv buffer.
    type(shr_horiz_remap_t), intent(in) :: rd
    integer, intent(in) :: gid
    integer :: idx
    idx = bsearch(rd%src_need_gids, rd%n_src_need, gid)
    src_gid_to_recvidx = rd%src_need_recvidx(idx)
  end function src_gid_to_recvidx

  !-------------------------------------------------------------------------------------------
  subroutine shr_horiz_remap_read_mapfile(rd, mapfile, comm, myrank, nprocs, &
       iosystem, ierr)
    !--------------------------------------------------------------------------
    ! Phase 1: Read the ESMF/TempestRemap mapping file, partition the target
    ! grid across MPI ranks, and extract the local sparse matrix entries.
    !
    ! After this call, rd%src_need_gids contains the sorted list of global
    ! source column IDs needed by this rank. The component wrapper must then
    ! fill gcol_to_rank(1:n_a) and call shr_horiz_remap_build_comm.
    !--------------------------------------------------------------------------
    use pio, only: file_desc_t, pio_closefile, pio_openfile, &
                   pio_nowrite, pio_inq_dimid, pio_inq_dimlen, &
                   pio_inq_varid, pio_get_var, var_desc_t, &
                   pio_noerr, PIO_IOTYPE_NETCDF, iosystem_desc_t

    type(shr_horiz_remap_t), intent(inout) :: rd
    character(len=*), intent(in) :: mapfile
    integer, intent(in) :: comm, myrank, nprocs
    type(iosystem_desc_t), intent(inout) :: iosystem
    integer, intent(out) :: ierr

    type(file_desc_t) :: pioid
    type(var_desc_t)  :: vid
    integer           :: dimid
    integer           :: n_s, n_a, n_b, dst_grid_rank
    integer, allocatable :: dst_grid_dims(:)
    integer, allocatable :: row_all(:), col_all(:)
    real(r8), allocatable :: S_all(:)
    real(r8), allocatable :: xc_b(:), yc_b(:)
    integer :: i, cnt, n_unique
    integer :: row_start, row_end, nlat, nlon
    integer, allocatable :: gcol_marker(:)

    ierr = 0

    ! Open map file via PIO
    ierr = pio_openfile(iosystem, pioid, PIO_IOTYPE_NETCDF, trim(mapfile), pio_nowrite)
    if (ierr /= pio_noerr) then
      ierr = 1
      return
    end if

    ! Read dimensions
    ierr = pio_inq_dimid(pioid, 'n_s', dimid)
    ierr = pio_inq_dimlen(pioid, dimid, n_s)
    ierr = pio_inq_dimid(pioid, 'n_a', dimid)
    ierr = pio_inq_dimlen(pioid, dimid, n_a)
    ierr = pio_inq_dimid(pioid, 'n_b', dimid)
    ierr = pio_inq_dimlen(pioid, dimid, n_b)

    rd%n_a = n_a
    rd%n_b = n_b

    ! Read dst_grid_rank and dims to determine nlat/nlon
    ierr = pio_inq_dimid(pioid, 'dst_grid_rank', dimid)
    ierr = pio_inq_dimlen(pioid, dimid, dst_grid_rank)
    if (dst_grid_rank /= 2) then
      ierr = 2
      call pio_closefile(pioid)
      return
    end if
    allocate(dst_grid_dims(dst_grid_rank))
    ierr = pio_inq_varid(pioid, 'dst_grid_dims', vid)
    ierr = pio_get_var(pioid, vid, dst_grid_dims)
    nlon = dst_grid_dims(1)
    nlat = dst_grid_dims(2)
    deallocate(dst_grid_dims)

    if (nlat * nlon /= n_b) then
      ierr = 3
      call pio_closefile(pioid)
      return
    end if
    rd%nlat = nlat
    rd%nlon = nlon

    ! Read full sparse matrix (all ranks read all - simple approach)
    allocate(row_all(n_s), col_all(n_s), S_all(n_s))
    ierr = pio_inq_varid(pioid, 'row', vid)
    ierr = pio_get_var(pioid, vid, row_all)
    ierr = pio_inq_varid(pioid, 'col', vid)
    ierr = pio_get_var(pioid, vid, col_all)
    ierr = pio_inq_varid(pioid, 'S', vid)
    ierr = pio_get_var(pioid, vid, S_all)

    ! Read target grid coordinates
    allocate(xc_b(n_b), yc_b(n_b))
    ierr = pio_inq_varid(pioid, 'xc_b', vid)
    ierr = pio_get_var(pioid, vid, xc_b)
    ierr = pio_inq_varid(pioid, 'yc_b', vid)
    ierr = pio_get_var(pioid, vid, yc_b)

    call pio_closefile(pioid)

    ! Extract unique lat/lon arrays from target coordinates
    ! Points are lon-major: point(i) has ilon = mod(i-1,nlon)+1, ilat = (i-1)/nlon+1
    allocate(rd%lon(nlon), rd%lat(nlat))
    do i = 1, nlon
      rd%lon(i) = xc_b(i)
    end do
    do i = 1, nlat
      rd%lat(i) = yc_b((i-1)*nlon + 1)
    end do
    deallocate(xc_b, yc_b)

    ! Partition target grid rows across MPI ranks (contiguous blocks)
    row_start = myrank * n_b / nprocs + 1
    row_end   = (myrank + 1) * n_b / nprocs
    rd%n_b_local = max(0, row_end - row_start + 1)
    rd%row_start = row_start

    ! Extract local sparse matrix entries (rows belonging to this rank)
    cnt = 0
    do i = 1, n_s
      if (row_all(i) >= row_start .and. row_all(i) <= row_end) cnt = cnt + 1
    end do
    rd%nnz_local = cnt

    allocate(rd%dst_local(cnt), rd%src_gid(cnt), rd%wgt(cnt))
    cnt = 0
    do i = 1, n_s
      if (row_all(i) >= row_start .and. row_all(i) <= row_end) then
        cnt = cnt + 1
        rd%dst_local(cnt) = row_all(i) - row_start + 1
        rd%src_gid(cnt) = col_all(i)
        rd%wgt(cnt) = S_all(i)
      end if
    end do
    deallocate(row_all, col_all, S_all)

    ! Find unique source columns needed by this rank
    allocate(gcol_marker(n_a))
    gcol_marker(:) = 0
    do i = 1, rd%nnz_local
      gcol_marker(rd%src_gid(i)) = 1
    end do

    n_unique = 0
    do i = 1, n_a
      if (gcol_marker(i) > 0) n_unique = n_unique + 1
    end do
    rd%n_src_need = n_unique

    ! Build sorted list of unique source column IDs
    allocate(rd%src_need_gids(n_unique), rd%src_need_recvidx(n_unique))
    cnt = 0
    do i = 1, n_a
      if (gcol_marker(i) > 0) then
        cnt = cnt + 1
        rd%src_need_gids(cnt) = i
      end if
    end do
    deallocate(gcol_marker)

    ierr = 0

  end subroutine shr_horiz_remap_read_mapfile

  !-------------------------------------------------------------------------------------------
  subroutine shr_horiz_remap_build_comm(rd, gcol_to_rank, comm, myrank, nprocs, &
       send_gcol_list, ierr)
    !--------------------------------------------------------------------------
    ! Phase 2: Build MPI communication pattern from the gcol_to_rank ownership
    ! array provided by the component.
    !
    ! Input:
    !   gcol_to_rank(1:n_a) - owning MPI rank for each global source column
    !                         (filled by component using its own decomposition)
    ! Output:
    !   send_gcol_list(1:n_send_total) - global column IDs that this rank must
    !                                     send. The component converts these to
    !                                     local indices for send buffer packing.
    !--------------------------------------------------------------------------
    use mpi, only: MPI_INTEGER

    type(shr_horiz_remap_t), intent(inout) :: rd
    integer, intent(in) :: gcol_to_rank(:)  ! (n_a)
    integer, intent(in) :: comm, myrank, nprocs
    integer, allocatable, intent(out) :: send_gcol_list(:)  ! (n_send_total)
    integer, intent(out) :: ierr

    integer :: i, r, gcol, owner_rank, cnt
    integer, allocatable :: need_from_rank(:), recv_gcols(:)

    ierr = 0

    ! Determine recv_counts: how many columns I need from each rank
    allocate(need_from_rank(0:nprocs-1))
    need_from_rank = 0
    do i = 1, rd%n_src_need
      gcol = rd%src_need_gids(i)
      owner_rank = gcol_to_rank(gcol)
      if (owner_rank >= 0 .and. owner_rank < nprocs) then
        need_from_rank(owner_rank) = need_from_rank(owner_rank) + 1
      end if
    end do

    allocate(rd%recv_counts(0:nprocs-1), rd%recv_displs(0:nprocs-1))
    rd%recv_counts = need_from_rank

    rd%recv_displs(0) = 0
    do r = 1, nprocs-1
      rd%recv_displs(r) = rd%recv_displs(r-1) + rd%recv_counts(r-1)
    end do
    rd%n_recv_total = sum(rd%recv_counts)

    ! Build ordered list of gcols to receive (grouped by source rank)
    allocate(recv_gcols(rd%n_recv_total))
    need_from_rank = rd%recv_displs  ! reuse as running offset
    do i = 1, rd%n_src_need
      gcol = rd%src_need_gids(i)
      owner_rank = gcol_to_rank(gcol)
      if (owner_rank >= 0 .and. owner_rank < nprocs) then
        need_from_rank(owner_rank) = need_from_rank(owner_rank) + 1
        recv_gcols(need_from_rank(owner_rank)) = gcol
      end if
    end do

    ! Build src_need_recvidx: for each src_need_gids entry, store its position
    ! in the recv buffer
    do i = 1, rd%n_src_need
      gcol = rd%src_need_gids(i)
      do cnt = 1, rd%n_recv_total
        if (recv_gcols(cnt) == gcol) then
          rd%src_need_recvidx(i) = cnt
          exit
        end if
      end do
    end do
    deallocate(need_from_rank)

    ! Exchange recv_counts so each rank knows what to send
    allocate(rd%send_counts(0:nprocs-1), rd%send_displs(0:nprocs-1))

    call mpi_alltoall(rd%recv_counts, 1, MPI_INTEGER, &
                      rd%send_counts, 1, MPI_INTEGER, &
                      comm, ierr)

    rd%send_displs(0) = 0
    do r = 1, nprocs-1
      rd%send_displs(r) = rd%send_displs(r-1) + rd%send_counts(r-1)
    end do
    rd%n_send_total = sum(rd%send_counts)

    ! Exchange the actual gcol IDs so senders know what to send
    allocate(send_gcol_list(rd%n_send_total))

    call mpi_alltoallv(recv_gcols, rd%recv_counts, rd%recv_displs, MPI_INTEGER, &
                       send_gcol_list, rd%send_counts, rd%send_displs, MPI_INTEGER, &
                       comm, ierr)

    deallocate(recv_gcols)

    rd%initialized = .true.

  end subroutine shr_horiz_remap_build_comm

  !-------------------------------------------------------------------------------------------
  subroutine shr_horiz_remap_apply(rd, send_buf, numlev, fld_out, comm, nprocs, ierr)
    !--------------------------------------------------------------------------
    ! Phase 3 (runtime): Apply horizontal remapping to a pre-packed send buffer.
    !
    ! The component packs send_buf(n_send_total * numlev) from its own data
    ! layout using the local indices it built from send_gcol_list. This routine
    ! does the MPI exchange and sparse matrix multiply.
    !
    ! Output: fld_out(n_b_local, numlev) - remapped field on lat-lon grid
    !--------------------------------------------------------------------------
    use mpi, only: MPI_DOUBLE_PRECISION

    type(shr_horiz_remap_t), intent(in) :: rd
    real(r8), intent(in)    :: send_buf(:)  ! (n_send_total * numlev)
    integer,  intent(in)    :: numlev
    real(r8), intent(out)   :: fld_out(:,:) ! (n_b_local, numlev) - preallocated by caller
    integer,  intent(in)    :: comm, nprocs
    integer,  intent(out)   :: ierr

    real(r8), allocatable :: recv_buf(:)
    integer, allocatable :: send_counts_lev(:), send_displs_lev(:)
    integer, allocatable :: recv_counts_lev(:), recv_displs_lev(:)
    integer :: i, k, src_local, dst_local

    ierr = 0
    fld_out = 0.0_r8

    allocate(recv_buf(rd%n_recv_total * numlev))

    ! Scale counts/displacements by numlev
    allocate(send_counts_lev(0:nprocs-1), send_displs_lev(0:nprocs-1))
    allocate(recv_counts_lev(0:nprocs-1), recv_displs_lev(0:nprocs-1))
    do i = 0, nprocs-1
      send_counts_lev(i) = rd%send_counts(i) * numlev
      send_displs_lev(i) = rd%send_displs(i) * numlev
      recv_counts_lev(i) = rd%recv_counts(i) * numlev
      recv_displs_lev(i) = rd%recv_displs(i) * numlev
    end do

    call mpi_alltoallv(send_buf, send_counts_lev, send_displs_lev, MPI_DOUBLE_PRECISION, &
                       recv_buf, recv_counts_lev, recv_displs_lev, MPI_DOUBLE_PRECISION, &
                       comm, ierr)

    deallocate(send_counts_lev, send_displs_lev, recv_counts_lev, recv_displs_lev)

    ! Sparse matrix-vector multiply
    do i = 1, rd%nnz_local
      dst_local = rd%dst_local(i)
      src_local = src_gid_to_recvidx(rd, rd%src_gid(i))
      do k = 1, numlev
        fld_out(dst_local, k) = fld_out(dst_local, k) + &
             rd%wgt(i) * recv_buf((src_local-1)*numlev + k)
      end do
    end do

    deallocate(recv_buf)

  end subroutine shr_horiz_remap_apply

end module shr_horiz_remap_mod
