module shr_horiz_remap_mod
  !-------------------------------------------------------------------------------------------
  !
  ! Shared horizontal remapping infrastructure for E3SM Fortran components.
  !
  ! Reads ESMF/TempestRemap mapping files and applies sparse matrix-vector
  ! multiply to remap fields from a source (neNpgP) grid to a target
  ! lat-lon grid.  Component-independent; each component provides a thin
  ! wrapper that handles its own decomposition and I/O.
  !
  ! Key design choices:
  !   - CRS sparse matrix (row_offsets, col_recvidx, weights)
  !   - Source-gather communication (alltoallv gathers needed source
  !     columns, then local CRS SpMV produces output)
  !   - Persistent MPI buffers (no per-call heap allocation)
  !   - First-contribution SpMV initialization (avoids zeroing output)
  !
  ! Usage (three phases):
  !   Phase 1: shr_horiz_remap_read_mapfile  - read map, partition target grid
  !   Phase 2: shr_horiz_remap_build_comm    - build MPI pattern + CRS matrix
  !   Phase 3: shr_horiz_remap_apply         - alltoallv + CRS SpMV
  !
  ! Component wrapper responsibilities:
  !   1. Calls phase 1 (shared: reads map file)
  !   2. Fills gcol_to_rank(1:n_a) using its own decomposition info
  !   3. Calls phase 2 (shared: builds MPI pattern + CRS, returns send_gcol_list)
  !   4. Converts send_gcol_list to local indices using its own indexing
  !   5. At runtime: packs send_buf from its own data layout, calls phase 3
  !
  !-------------------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private
  save

  public :: shr_horiz_remap_t
  public :: SHR_FILL_VALUE

  ! Standard fill value for remapped output.  Matches CAM's fillvalue
  ! (cam_history_support) and is representable in both real(4) and real(8).
  real(r8), parameter :: SHR_FILL_VALUE = 1.0e20_r8

  type shr_horiz_remap_t
    logical  :: initialized = .false.
    integer  :: n_a = 0          ! source grid size (global)
    integer  :: n_b = 0          ! target grid size (global)
    integer  :: n_b_local = 0    ! target points owned by this rank
    integer  :: row_start = 0    ! global starting target row (1-indexed)

    ! Target lat-lon grid dimensions and coordinates
    integer  :: nlat = 0
    integer  :: nlon = 0
    real(r8), allocatable :: lat(:)   ! (nlat)
    real(r8), allocatable :: lon(:)   ! (nlon)

    ! Target cell ocean fraction for renormalization of remapped values.
    ! frac_b_local(i) = sum of weights for target cell i; dividing by this
    ! corrects for partial ocean coverage at coastlines.
    real(r8), allocatable :: frac_b_local(:)  ! (n_b_local)

    ! CRS sparse matrix for SpMV (built during build_comm).
    ! All entries for a given output row are contiguous — cache-friendly
    ! for the SpMV inner loop.
    integer  :: nnz_local = 0
    integer,  allocatable :: row_offsets(:)   ! (n_b_local+1) CRS row pointers (1-indexed)
    integer,  allocatable :: col_recvidx(:)   ! (nnz_local) recv-buffer position per nonzero
    real(r8), allocatable :: wgt(:)           ! (nnz_local) interpolation weights

    ! Temporary triplet data: populated by read_mapfile, consumed by build_comm
    integer,  allocatable :: tmp_dst(:)       ! local dest row (1..n_b_local)
    integer,  allocatable :: tmp_src_gid(:)   ! global source column ID
    real(r8), allocatable :: tmp_wgt(:)       ! weight
    integer  :: tmp_nnz = 0

    ! Source column gathering (init-only; freed at end of build_comm)
    integer  :: n_src_need = 0
    integer,  allocatable :: src_need_gids(:) ! sorted unique global IDs of needed src cols

    ! MPI communication pattern for Alltoallv (source-gather)
    integer,  allocatable :: send_counts(:)   ! (0:npes-1)
    integer,  allocatable :: send_displs(:)
    integer,  allocatable :: recv_counts(:)
    integer,  allocatable :: recv_displs(:)
    integer  :: n_send_total = 0
    integer  :: n_recv_total = 0

    ! Persistent workspace for apply (grows as needed)
    real(r8), allocatable :: ws_recv_buf(:)
  contains
    procedure :: read_mapfile => shr_horiz_remap_read_mapfile
    procedure :: build_comm   => shr_horiz_remap_build_comm
    procedure :: apply        => shr_horiz_remap_apply
    procedure :: apply_masked => shr_horiz_remap_apply_masked
  end type shr_horiz_remap_t

CONTAINS

  !-------------------------------------------------------------------------------------------
  subroutine shr_horiz_remap_read_mapfile(rd, mapfile, comm, myrank, nprocs, &
       iosystem, ierr)
    !--------------------------------------------------------------------------
    ! Phase 1: Read the ESMF/TempestRemap mapping file, partition the target
    ! grid across MPI ranks, and extract local triplets into temporary arrays.
    ! Call build_comm next.
    !--------------------------------------------------------------------------
    use pio, only: file_desc_t, pio_closefile, pio_openfile, &
                   pio_nowrite, pio_inq_dimid, pio_inq_dimlen, &
                   pio_inq_varid, pio_get_var, var_desc_t, &
                   pio_noerr, PIO_IOTYPE_NETCDF, iosystem_desc_t

    class(shr_horiz_remap_t), intent(inout) :: rd
    character(len=*), intent(in) :: mapfile
    integer, intent(in) :: comm, myrank, nprocs
    type(iosystem_desc_t), intent(inout) :: iosystem
    integer, intent(out) :: ierr

    type(file_desc_t) :: pioid
    integer :: i, cnt, n_unique
    integer :: n_s, n_a, n_b, row_start, row_end
    integer, allocatable :: row_all(:), col_all(:)
    real(r8), allocatable :: S_all(:)
    integer, allocatable :: gcol_marker(:)

    ierr = 0

    ! Guard against double-initialization (e.g. compute_on_startup calling init twice)
    if (rd%initialized) return

    ierr = pio_openfile(iosystem, pioid, PIO_IOTYPE_NETCDF, trim(mapfile), pio_nowrite)
    if (ierr /= pio_noerr) then
      ierr = 1; return
    end if
    call read_pio_contents(pioid, n_s, n_a, n_b, row_all, col_all, S_all, ierr)
    if (ierr /= 0) return

    rd%n_a = n_a
    rd%n_b = n_b

    ! Partition target grid across ranks (contiguous blocks)
    row_start = myrank * n_b / nprocs + 1
    row_end   = (myrank + 1) * n_b / nprocs
    rd%n_b_local = max(0, row_end - row_start + 1)
    rd%row_start = row_start

    ! Extract local triplets (target rows belonging to this rank)
    cnt = 0
    do i = 1, n_s
      if (row_all(i) >= row_start .and. row_all(i) <= row_end) cnt = cnt + 1
    end do
    rd%tmp_nnz = cnt

    allocate(rd%tmp_dst(cnt), rd%tmp_src_gid(cnt), rd%tmp_wgt(cnt))
    cnt = 0
    do i = 1, n_s
      if (row_all(i) >= row_start .and. row_all(i) <= row_end) then
        cnt = cnt + 1
        rd%tmp_dst(cnt)     = row_all(i) - row_start + 1
        rd%tmp_src_gid(cnt) = col_all(i)
        rd%tmp_wgt(cnt)     = S_all(i)
      end if
    end do
    deallocate(row_all, col_all, S_all)

    ! Build sorted list of unique source columns needed by this rank
    allocate(gcol_marker(n_a))
    gcol_marker(:) = 0
    do i = 1, rd%tmp_nnz
      gcol_marker(rd%tmp_src_gid(i)) = 1
    end do

    n_unique = 0
    do i = 1, n_a
      if (gcol_marker(i) > 0) n_unique = n_unique + 1
    end do
    rd%n_src_need = n_unique

    allocate(rd%src_need_gids(n_unique))
    cnt = 0
    do i = 1, n_a
      if (gcol_marker(i) > 0) then
        cnt = cnt + 1
        rd%src_need_gids(cnt) = i
      end if
    end do
    deallocate(gcol_marker)

  contains

    subroutine read_pio_contents(pioid, n_s, n_a, n_b, row_all, col_all, S_all, ierr)
      !------------------------------------------------------------------------
      ! Read all required data from the map file.
      ! Always closes pioid before returning.
      !------------------------------------------------------------------------
      type(file_desc_t), intent(inout) :: pioid
      integer, intent(out) :: n_s, n_a, n_b, ierr
      integer, allocatable, intent(out) :: row_all(:), col_all(:)
      real(r8), allocatable, intent(out) :: S_all(:)

      type(var_desc_t) :: vid
      integer :: dimid, pio_ierr, dst_grid_rank, nlat, nlon, i
      integer :: dst_grid_dims(2)
      real(r8), allocatable :: xc_b(:), yc_b(:), frac_b_all(:)

      ierr = 0

      ! --- dimensions ---
      pio_ierr = pio_inq_dimid(pioid, 'n_s', dimid)
      if (pio_ierr /= pio_noerr) then; ierr = 10; call pio_closefile(pioid); return; end if
      pio_ierr = pio_inq_dimlen(pioid, dimid, n_s)

      pio_ierr = pio_inq_dimid(pioid, 'n_a', dimid)
      if (pio_ierr /= pio_noerr) then; ierr = 11; call pio_closefile(pioid); return; end if
      pio_ierr = pio_inq_dimlen(pioid, dimid, n_a)

      pio_ierr = pio_inq_dimid(pioid, 'n_b', dimid)
      if (pio_ierr /= pio_noerr) then; ierr = 12; call pio_closefile(pioid); return; end if
      pio_ierr = pio_inq_dimlen(pioid, dimid, n_b)

      ! --- target grid must be structured lat-lon ---
      pio_ierr = pio_inq_dimid(pioid, 'dst_grid_rank', dimid)
      if (pio_ierr /= pio_noerr) then; ierr = 13; call pio_closefile(pioid); return; end if
      pio_ierr = pio_inq_dimlen(pioid, dimid, dst_grid_rank)
      if (dst_grid_rank /= 2) then; ierr = 2; call pio_closefile(pioid); return; end if

      pio_ierr = pio_inq_varid(pioid, 'dst_grid_dims', vid)
      if (pio_ierr /= pio_noerr) then; ierr = 14; call pio_closefile(pioid); return; end if
      pio_ierr = pio_get_var(pioid, vid, dst_grid_dims)
      nlon = dst_grid_dims(1); nlat = dst_grid_dims(2)
      if (nlat * nlon /= n_b) then; ierr = 3; call pio_closefile(pioid); return; end if
      rd%nlat = nlat; rd%nlon = nlon

      ! --- sparse matrix triplets ---
      allocate(row_all(n_s), col_all(n_s), S_all(n_s))
      pio_ierr = pio_inq_varid(pioid, 'row', vid)
      if (pio_ierr /= pio_noerr) then; ierr = 15; call pio_closefile(pioid); return; end if
      pio_ierr = pio_get_var(pioid, vid, row_all)
      pio_ierr = pio_inq_varid(pioid, 'col', vid)
      if (pio_ierr /= pio_noerr) then; ierr = 16; call pio_closefile(pioid); return; end if
      pio_ierr = pio_get_var(pioid, vid, col_all)
      pio_ierr = pio_inq_varid(pioid, 'S', vid)
      if (pio_ierr /= pio_noerr) then; ierr = 17; call pio_closefile(pioid); return; end if
      pio_ierr = pio_get_var(pioid, vid, S_all)

      ! --- target grid coordinates ---
      allocate(xc_b(n_b), yc_b(n_b))
      pio_ierr = pio_inq_varid(pioid, 'xc_b', vid)
      if (pio_ierr /= pio_noerr) then; ierr = 18; call pio_closefile(pioid); return; end if
      pio_ierr = pio_get_var(pioid, vid, xc_b)
      pio_ierr = pio_inq_varid(pioid, 'yc_b', vid)
      if (pio_ierr /= pio_noerr) then; ierr = 19; call pio_closefile(pioid); return; end if
      pio_ierr = pio_get_var(pioid, vid, yc_b)

      ! frac_b (target cell coverage fraction) is computed in build_comm
      ! from the CRS weight sum directly, avoiding PIO state corruption.

      call pio_closefile(pioid)

      ! Extract 1D lon/lat axes (lon-major ordering)
      allocate(rd%lon(nlon), rd%lat(nlat))
      rd%lon(1:nlon) = xc_b(1:nlon)
      do i = 1, nlat
        rd%lat(i) = yc_b((i-1)*nlon + 1)
      end do
      deallocate(xc_b, yc_b)

    end subroutine read_pio_contents

  end subroutine shr_horiz_remap_read_mapfile

  !-------------------------------------------------------------------------------------------
  subroutine shr_horiz_remap_build_comm(rd, gcol_to_rank, comm, myrank, nprocs, &
       send_gcol_list, ierr)
    !--------------------------------------------------------------------------
    ! Phase 2: Build MPI communication pattern and convert temporary triplets
    ! into CRS format for the SpMV.
    !
    ! Source-gather: each rank tells other ranks which source columns it needs.
    ! After alltoallv, each rank has the source data for its local SpMV.
    !
    ! Input:
    !   gcol_to_rank(1:n_a) - owning MPI rank for each global source column
    ! Output:
    !   send_gcol_list(1:n_send_total) - global column IDs this rank must send
    !--------------------------------------------------------------------------
    use mpi, only: MPI_INTEGER, MPI_SUCCESS

    class(shr_horiz_remap_t), intent(inout) :: rd
    integer, intent(in) :: gcol_to_rank(:)  ! (n_a)
    integer, intent(in) :: comm, myrank, nprocs
    integer, allocatable, intent(out) :: send_gcol_list(:)
    integer, intent(out) :: ierr

    integer :: i, r, gcol, owner_rank, irow, j
    integer, allocatable :: need_from_rank(:), recv_gcols(:)
    integer, allocatable :: gcol_to_recvpos(:)
    integer, allocatable :: recvidx_unsorted(:), row_counts(:), bucket_pos(:)

    ierr = 0

    ! Guard against double-initialization (e.g. compute_on_startup calling init twice)
    if (rd%initialized) return

    ! --- Step 1: determine recv_counts ---
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

    ! --- Step 2: build ordered list of gcols to receive (grouped by source rank) ---
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
    deallocate(need_from_rank)

    ! --- Step 3: map each unique source gcol to its recv buffer position ---
    ! Use a temporary array indexed by global column ID: O(n_a) memory, O(n_recv) time
    allocate(gcol_to_recvpos(rd%n_a))
    gcol_to_recvpos(:) = 0
    do i = 1, rd%n_recv_total
      gcol_to_recvpos(recv_gcols(i)) = i
    end do

    ! --- Step 4: exchange counts and gcol IDs ---
    allocate(rd%send_counts(0:nprocs-1), rd%send_displs(0:nprocs-1))

    call mpi_alltoall(rd%recv_counts, 1, MPI_INTEGER, &
                      rd%send_counts, 1, MPI_INTEGER, comm, ierr)
    if (ierr /= MPI_SUCCESS) return

    rd%send_displs(0) = 0
    do r = 1, nprocs-1
      rd%send_displs(r) = rd%send_displs(r-1) + rd%send_counts(r-1)
    end do
    rd%n_send_total = sum(rd%send_counts)

    allocate(send_gcol_list(rd%n_send_total))
    call mpi_alltoallv(recv_gcols, rd%recv_counts, rd%recv_displs, MPI_INTEGER, &
                       send_gcol_list, rd%send_counts, rd%send_displs, MPI_INTEGER, &
                       comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    deallocate(recv_gcols)

    ! --- Step 5: build CRS matrix from temporary triplets ---
    rd%nnz_local = rd%tmp_nnz

    ! Map each nonzero's source gcol to its recv buffer position
    allocate(recvidx_unsorted(rd%nnz_local))
    do i = 1, rd%nnz_local
      recvidx_unsorted(i) = gcol_to_recvpos(rd%tmp_src_gid(i))
    end do
    deallocate(gcol_to_recvpos)

    ! Bucket-sort by destination row to build CRS row_offsets
    allocate(row_counts(rd%n_b_local))
    row_counts = 0
    do i = 1, rd%nnz_local
      row_counts(rd%tmp_dst(i)) = row_counts(rd%tmp_dst(i)) + 1
    end do

    allocate(rd%row_offsets(rd%n_b_local + 1))
    rd%row_offsets(1) = 1
    do i = 1, rd%n_b_local
      rd%row_offsets(i+1) = rd%row_offsets(i) + row_counts(i)
    end do

    ! Scatter entries into CRS order
    allocate(rd%col_recvidx(rd%nnz_local), rd%wgt(rd%nnz_local))
    allocate(bucket_pos(rd%n_b_local))
    bucket_pos(1:rd%n_b_local) = rd%row_offsets(1:rd%n_b_local)

    do i = 1, rd%nnz_local
      irow = rd%tmp_dst(i)
      j = bucket_pos(irow)
      rd%col_recvidx(j) = recvidx_unsorted(i)
      rd%wgt(j) = rd%tmp_wgt(i)
      bucket_pos(irow) = bucket_pos(irow) + 1
    end do

    deallocate(recvidx_unsorted, row_counts, bucket_pos)

    ! Compute target cell coverage fraction for coastal renormalization.
    ! frac_b = sum of remap weights per target row. When the source grid
    ! excludes land cells, frac_b < 1 at coastlines; dividing the SpMV
    ! result by frac_b corrects for partial coverage.
    allocate(rd%frac_b_local(rd%n_b_local))
    rd%frac_b_local(:) = 0.0_r8
    do i = 1, rd%n_b_local
      do j = rd%row_offsets(i), rd%row_offsets(i+1) - 1
        rd%frac_b_local(i) = rd%frac_b_local(i) + rd%wgt(j)
      end do
    end do

    ! Free init-only data
    deallocate(rd%tmp_dst, rd%tmp_src_gid, rd%tmp_wgt)
    rd%tmp_nnz = 0
    deallocate(rd%src_need_gids)
    rd%n_src_need = 0

    rd%initialized = .true.

  end subroutine shr_horiz_remap_build_comm

  !-------------------------------------------------------------------------------------------
  subroutine shr_horiz_remap_apply(rd, send_buf, numlev, fld_out, comm, nprocs, ierr)
    !--------------------------------------------------------------------------
    ! Phase 3 (runtime): Source-gather remap.
    !   1. Alltoallv gathers needed source column data
    !   2. CRS SpMV produces output on local target rows
    !
    ! Input:  send_buf(n_send_total * numlev) - packed source data
    ! Output: fld_out(n_b_local, numlev)
    !--------------------------------------------------------------------------
    use mpi, only: MPI_DOUBLE_PRECISION

    class(shr_horiz_remap_t), intent(inout) :: rd
    real(r8), intent(in)    :: send_buf(:)
    integer,  intent(in)    :: numlev
    real(r8), intent(out)   :: fld_out(:,:)
    integer,  intent(in)    :: comm, nprocs
    integer,  intent(out)   :: ierr

    integer :: needed, i, k, j, src_idx, jbeg, jend
    integer :: send_counts_lev(0:nprocs-1), send_displs_lev(0:nprocs-1)
    integer :: recv_counts_lev(0:nprocs-1), recv_displs_lev(0:nprocs-1)

    ierr = 0

    ! Grow persistent recv workspace if needed.
    ! Always allocate at least size 1 so MPI_Alltoallv never receives
    ! an unallocated array on ranks with no target cells.
    needed = rd%n_recv_total * numlev
    if (.not. allocated(rd%ws_recv_buf) .or. size(rd%ws_recv_buf) < max(1, needed)) then
      if (allocated(rd%ws_recv_buf)) deallocate(rd%ws_recv_buf)
      allocate(rd%ws_recv_buf(max(1, needed)))
    end if

    ! Scale counts/displacements by numlev
    do i = 0, nprocs-1
      send_counts_lev(i) = rd%send_counts(i) * numlev
      send_displs_lev(i) = rd%send_displs(i) * numlev
      recv_counts_lev(i) = rd%recv_counts(i) * numlev
      recv_displs_lev(i) = rd%recv_displs(i) * numlev
    end do

    call mpi_alltoallv(send_buf, send_counts_lev, send_displs_lev, MPI_DOUBLE_PRECISION, &
                       rd%ws_recv_buf, recv_counts_lev, recv_displs_lev, MPI_DOUBLE_PRECISION, &
                       comm, ierr)
    if (ierr /= 0) return

    ! CRS SpMV with first-contribution initialization
    do i = 1, rd%n_b_local
      jbeg = rd%row_offsets(i)
      jend = rd%row_offsets(i+1) - 1
      if (jbeg <= jend) then
        ! First contribution: initialize
        src_idx = rd%col_recvidx(jbeg)
        do k = 1, numlev
          fld_out(i, k) = rd%wgt(jbeg) * rd%ws_recv_buf((src_idx-1)*numlev + k)
        end do
        ! Remaining: accumulate
        do j = jbeg+1, jend
          src_idx = rd%col_recvidx(j)
          do k = 1, numlev
            fld_out(i, k) = fld_out(i, k) + &
                 rd%wgt(j) * rd%ws_recv_buf((src_idx-1)*numlev + k)
          end do
        end do
      else
        do k = 1, numlev
          fld_out(i, k) = 0.0_r8
        end do
      end if
    end do

    ! Renormalize by frac_b to correct for partial ocean coverage.
    ! Without this, coastal target cells are systematically low because
    ! the weight sum < 1 (land source cells contribute zero).
    if (allocated(rd%frac_b_local)) then
      do i = 1, rd%n_b_local
        if (rd%frac_b_local(i) > 1.0e-6_r8) then
          do k = 1, numlev
            fld_out(i, k) = fld_out(i, k) / rd%frac_b_local(i)
          end do
        else
          ! Pure land cell: mark as fill
          do k = 1, numlev
            fld_out(i, k) = SHR_FILL_VALUE
          end do
        end if
      end do
    end if

  end subroutine shr_horiz_remap_apply

  !-------------------------------------------------------------------------------------------
  subroutine shr_horiz_remap_apply_masked(rd, send_buf, numlev, fld_out, comm, nprocs, ierr)
    !--------------------------------------------------------------------------
    ! Phase 3 variant: Source-gather remap with validity mask.
    !
    ! The caller packs send_buf with (numlev+1) entries per source cell:
    !   entries 1..numlev  = field values (fills zeroed to 0)
    !   entry   numlev+1   = validity mask (1.0 for valid, 0.0 for fill)
    !
    ! After SpMV, the mask level gives valid_frac(i) = sum(w_j * mask_j),
    ! the correct normalization denominator that excludes fill sources.
    ! This replaces frac_b normalization for properly masked output.
    !
    ! Input:  send_buf(n_send_total * (numlev+1))
    ! Output: fld_out(n_b_local, numlev) with SHR_FILL_VALUE where invalid
    !--------------------------------------------------------------------------
    use mpi, only: MPI_DOUBLE_PRECISION

    class(shr_horiz_remap_t), intent(inout) :: rd
    real(r8), intent(in)    :: send_buf(:)
    integer,  intent(in)    :: numlev
    real(r8), intent(out)   :: fld_out(:,:)
    integer,  intent(in)    :: comm, nprocs
    integer,  intent(out)   :: ierr

    integer :: nlev_packed, needed, i, k, j, src_idx, jbeg, jend
    integer :: send_counts_lev(0:nprocs-1), send_displs_lev(0:nprocs-1)
    integer :: recv_counts_lev(0:nprocs-1), recv_displs_lev(0:nprocs-1)
    real(r8) :: valid_frac
    real(r8), parameter :: vfrac_eps = 1.0e-10_r8

    ierr = 0
    nlev_packed = numlev + 1

    ! Grow persistent recv workspace if needed
    needed = rd%n_recv_total * nlev_packed
    if (.not. allocated(rd%ws_recv_buf) .or. size(rd%ws_recv_buf) < max(1, needed)) then
      if (allocated(rd%ws_recv_buf)) deallocate(rd%ws_recv_buf)
      allocate(rd%ws_recv_buf(max(1, needed)))
    end if

    ! Scale counts/displacements by nlev_packed
    do i = 0, nprocs-1
      send_counts_lev(i) = rd%send_counts(i) * nlev_packed
      send_displs_lev(i) = rd%send_displs(i) * nlev_packed
      recv_counts_lev(i) = rd%recv_counts(i) * nlev_packed
      recv_displs_lev(i) = rd%recv_displs(i) * nlev_packed
    end do

    call mpi_alltoallv(send_buf, send_counts_lev, send_displs_lev, MPI_DOUBLE_PRECISION, &
                       rd%ws_recv_buf, recv_counts_lev, recv_displs_lev, MPI_DOUBLE_PRECISION, &
                       comm, ierr)
    if (ierr /= 0) return

    ! CRS SpMV over field levels + mask level, with per-row normalization
    do i = 1, rd%n_b_local
      jbeg = rd%row_offsets(i)
      jend = rd%row_offsets(i+1) - 1
      if (jbeg <= jend) then
        ! First contribution
        src_idx = rd%col_recvidx(jbeg)
        do k = 1, numlev
          fld_out(i, k) = rd%wgt(jbeg) * rd%ws_recv_buf((src_idx-1)*nlev_packed + k)
        end do
        valid_frac = rd%wgt(jbeg) * rd%ws_recv_buf((src_idx-1)*nlev_packed + nlev_packed)
        ! Remaining contributions
        do j = jbeg+1, jend
          src_idx = rd%col_recvidx(j)
          do k = 1, numlev
            fld_out(i, k) = fld_out(i, k) + &
                 rd%wgt(j) * rd%ws_recv_buf((src_idx-1)*nlev_packed + k)
          end do
          valid_frac = valid_frac + &
               rd%wgt(j) * rd%ws_recv_buf((src_idx-1)*nlev_packed + nlev_packed)
        end do
        ! Normalize by valid_frac (the correct denominator excluding fills)
        if (valid_frac > vfrac_eps) then
          do k = 1, numlev
            fld_out(i, k) = fld_out(i, k) / valid_frac
          end do
        else
          do k = 1, numlev
            fld_out(i, k) = SHR_FILL_VALUE
          end do
        end if
      else
        ! No source contributions (target cell outside source domain)
        do k = 1, numlev
          fld_out(i, k) = SHR_FILL_VALUE
        end do
      end if
    end do

  end subroutine shr_horiz_remap_apply_masked

end module shr_horiz_remap_mod
