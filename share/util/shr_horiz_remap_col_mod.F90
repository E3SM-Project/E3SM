module shr_horiz_remap_col_mod
  !-------------------------------------------------------------------------------------------
  !
  ! Generic driver for shr_horiz_remap_mod for components whose horizontal
  ! decomposition is a plain list of global column indices.
  !
  ! EAM owns a bespoke wrapper (components/eam/src/control/horiz_remap_mod.F90)
  ! because its physics decomposition is (chunk, column-in-chunk) and its output
  ! goes through cam_grid / cam_history.  ELM and MOSART are much simpler: each
  ! rank owns a contiguous slice of a 1-D array whose entries carry a global
  ! index into the source grid (ldecomp%gdc2glo for ELM, rtmCTL%gindex for
  ! MOSART).  This module captures everything those two share:
  !
  !   * building gcol_to_rank from the caller's global-index list
  !   * converting the shared layer's send_gcol_list into local array offsets
  !   * packing the send buffer (field values + validity mask)
  !   * caching PIO decompositions and writing remapped fields
  !
  ! Usage:
  !   call remap%init(mapfile, my_gindex, n_a_expected, mpicom, iam, npes, &
  !                   iosystem, ierr, errmsg)
  !   call remap%remap_field(fld_local, numlev, fld_out)   ! fld_out allocated here
  !   call remap%write_field(File, varid, fld_out, numlev, PIO_DOUBLE)
  !
  ! Partial source coverage
  ! -----------------------
  ! ELM's land decomposition covers only the land cells of the source grid, yet
  ! the map file is built against the full (global) grid, so a target cell can
  ! need source columns that no rank owns.  shr_horiz_remap_build_comm refuses
  ! such a map outright (a dropped column would read out of bounds later), so
  ! this module assigns every unowned source column a deterministic *virtual*
  ! owner, mod(gcol-1, npes), spreading them evenly.  A virtual column has
  ! local index 0 and is packed according to `missing_as_zero`:
  !
  !   missing_as_zero = .false. (default)  value 0, mask 0 -> excluded from the
  !       average.  A target cell that is half land gets the mean over its land
  !       part, which is what a land emulator wants from QRUNOFF or TSOI.
  !   missing_as_zero = .true.             value 0, mask 1 -> counted as a real
  !       zero.  Remapping the land fraction this way yields the true land
  !       fraction of the target cell rather than 1.0 everywhere any land
  !       appears.
  !
  ! Fill values inside the caller's own data (spval columns) follow the same
  ! rule via the optional `fillval` argument.
  !
  !-------------------------------------------------------------------------------------------

  use shr_kind_mod,        only: r8 => shr_kind_r8
  use shr_horiz_remap_mod, only: shr_horiz_remap_t, SHR_FILL_VALUE
  use pio,                 only: io_desc_t, iosystem_desc_t

  implicit none
  private
  save

  public :: shr_horiz_remap_col_t
  public :: SHR_FILL_VALUE

  type shr_horiz_remap_col_t
    type(shr_horiz_remap_t) :: shared

    ! Local index (1-based, into the caller's packed local array) of each entry
    ! in the shared layer's send list.  Zero marks a virtual column -- one this
    ! rank was assigned only because no rank actually owns it.
    integer, allocatable :: send_local_idx(:)

    ! Persistent send workspace (grows as needed)
    real(r8), allocatable :: ws_send_buf(:)

    ! Cached PIO decompositions
    logical :: iodesc_2d_valid = .false.
    integer :: iodesc_2d_dtype = 0
    type(io_desc_t) :: iodesc_2d
    logical :: iodesc_3d_valid = .false.
    integer :: iodesc_3d_nlev  = 0
    integer :: iodesc_3d_dtype = 0
    type(io_desc_t) :: iodesc_3d

    type(iosystem_desc_t), pointer :: iosystem => null()
    integer :: mpicom = 0
    integer :: npes   = 0
  contains
    procedure :: init        => shr_horiz_remap_col_init
    procedure :: remap_field => shr_horiz_remap_col_field
    procedure :: write_field => shr_horiz_remap_col_write
    procedure :: is_active   => shr_horiz_remap_col_is_active
    procedure :: nlat        => shr_horiz_remap_col_nlat
    procedure :: nlon        => shr_horiz_remap_col_nlon
  end type shr_horiz_remap_col_t

CONTAINS

  !-------------------------------------------------------------------------------------------
  logical function shr_horiz_remap_col_is_active(self)
    class(shr_horiz_remap_col_t), intent(in) :: self
    shr_horiz_remap_col_is_active = self%shared%initialized
  end function shr_horiz_remap_col_is_active

  !-------------------------------------------------------------------------------------------
  integer function shr_horiz_remap_col_nlat(self)
    class(shr_horiz_remap_col_t), intent(in) :: self
    shr_horiz_remap_col_nlat = self%shared%nlat
  end function shr_horiz_remap_col_nlat

  !-------------------------------------------------------------------------------------------
  integer function shr_horiz_remap_col_nlon(self)
    class(shr_horiz_remap_col_t), intent(in) :: self
    shr_horiz_remap_col_nlon = self%shared%nlon
  end function shr_horiz_remap_col_nlon

  !-------------------------------------------------------------------------------------------
  subroutine shr_horiz_remap_col_init(self, mapfile, my_gindex, n_a_expected, &
       mpicom, iam, npes, iosystem, ierr, errmsg)
    !
    ! my_gindex(i) is the global source-grid index of the caller's i-th local
    ! column, using the same 1..n_a numbering as the map file's col_a.
    !
    use mpi, only: MPI_INTEGER, MPI_MAX

    class(shr_horiz_remap_col_t), intent(inout) :: self
    character(len=*),      intent(in)  :: mapfile
    integer,               intent(in)  :: my_gindex(:)
    integer,               intent(in)  :: n_a_expected
    integer,               intent(in)  :: mpicom, iam, npes
    type(iosystem_desc_t), pointer     :: iosystem
    integer,               intent(out) :: ierr
    character(len=*),      intent(out) :: errmsg

    integer :: i, g, n_local, mpierr
    integer, allocatable :: gcol_to_rank(:), my_owner(:)
    integer, allocatable :: gcol_to_local(:), send_gcol_list(:)

    ierr   = 0
    errmsg = ''

    self%iosystem => iosystem
    self%mpicom   = mpicom
    self%npes     = npes

    n_local = size(my_gindex)

    ! Phase 1: read the map file
    call self%shared%read_mapfile(mapfile, mpicom, iam, npes, iosystem, ierr)
    if (ierr /= 0) then
      errmsg = 'error reading mapping file '//trim(mapfile)
      return
    end if

    ! The masked-renormalization path is only well conditioned for non-negative
    ! maps; a negative weight lets valid_frac cancel at a partially covered cell
    ! and blows the result far outside the source range.  Refuse it here rather
    ! than shipping a corrupt tape.
    if (self%shared%has_negative_weights) then
      ierr = 1
      errmsg = 'mapping file '//trim(mapfile)//' has NEGATIVE weights, which '// &
           'are unsafe for the masked remap path; use a first-order '// &
           'conservative or bilinear map'
      return
    end if

    if (self%shared%n_a /= n_a_expected) then
      ierr = 1
      errmsg = 'mapping file '//trim(mapfile)//' n_a does not match the '// &
           'component source grid size'
      return
    end if

    ! Build gcol_to_rank.  Ownership is exclusive, so an MPI_MAX reduction over
    ! an array seeded with -1 picks out the single real owner of each column.
    allocate(my_owner(self%shared%n_a), gcol_to_rank(self%shared%n_a))
    my_owner(:) = -1
    do i = 1, n_local
      g = my_gindex(i)
      if (g < 1 .or. g > self%shared%n_a) then
        deallocate(my_owner, gcol_to_rank)
        ierr = 1
        errmsg = 'local global-index out of range for the source grid of '// &
             trim(mapfile)
        return
      end if
      my_owner(g) = iam
    end do

    call mpi_allreduce(my_owner, gcol_to_rank, self%shared%n_a, &
         MPI_INTEGER, MPI_MAX, mpicom, mpierr)
    deallocate(my_owner)
    if (mpierr /= 0) then
      deallocate(gcol_to_rank)
      ierr = 1
      errmsg = 'MPI_Allreduce failed while building the remap ownership map'
      return
    end if

    ! Assign a virtual owner to every column no rank actually holds (see the
    ! module header).  Round-robin keeps the extra send volume balanced.
    do g = 1, self%shared%n_a
      if (gcol_to_rank(g) < 0) gcol_to_rank(g) = mod(g - 1, npes)
    end do

    ! Phase 2: communication pattern + CRS matrix
    call self%shared%build_comm(gcol_to_rank, mpicom, iam, npes, send_gcol_list, ierr)
    deallocate(gcol_to_rank)
    if (ierr /= 0) then
      errmsg = 'error building the remap communication pattern for '//trim(mapfile)
      return
    end if

    ! Translate the send list into local array offsets.
    allocate(gcol_to_local(self%shared%n_a))
    gcol_to_local(:) = 0
    do i = 1, n_local
      gcol_to_local(my_gindex(i)) = i
    end do

    allocate(self%send_local_idx(max(1, self%shared%n_send_total)))
    self%send_local_idx(:) = 0
    do i = 1, self%shared%n_send_total
      self%send_local_idx(i) = gcol_to_local(send_gcol_list(i))
    end do

    deallocate(gcol_to_local)
    if (allocated(send_gcol_list)) deallocate(send_gcol_list)

  end subroutine shr_horiz_remap_col_init

  !-------------------------------------------------------------------------------------------
  subroutine shr_horiz_remap_col_field(self, fld_local, numlev, fld_out, &
       fillval, missing_as_zero)
    !
    ! Remap fld_local(n_local, numlev) onto the target grid.  fld_out is
    ! allocated here as (n_b_local, numlev) and carries SHR_FILL_VALUE wherever
    ! the target cell has no valid source coverage.
    !
    class(shr_horiz_remap_col_t), intent(inout) :: self
    real(r8),             intent(in)  :: fld_local(:,:)
    integer,              intent(in)  :: numlev
    real(r8), allocatable, intent(out) :: fld_out(:,:)
    real(r8), optional,   intent(in)  :: fillval
    logical,  optional,   intent(in)  :: missing_as_zero

    integer  :: i, k, idx, needed, nlev_packed, base, ierr
    real(r8) :: lfill, missing_mask
    logical  :: have_fill, is_missing

    allocate(fld_out(self%shared%n_b_local, numlev))
    fld_out(:,:) = SHR_FILL_VALUE

    if (.not. self%shared%initialized) return

    have_fill = present(fillval)
    lfill     = SHR_FILL_VALUE
    if (have_fill) lfill = fillval

    missing_mask = 0.0_r8
    if (present(missing_as_zero)) then
      if (missing_as_zero) missing_mask = 1.0_r8
    end if

    nlev_packed = numlev + 1   ! field levels + validity mask

    needed = max(1, self%shared%n_send_total * nlev_packed)
    if (.not. allocated(self%ws_send_buf) .or. size(self%ws_send_buf) < needed) then
      if (allocated(self%ws_send_buf)) deallocate(self%ws_send_buf)
      allocate(self%ws_send_buf(needed))
    end if

    do i = 1, self%shared%n_send_total
      base = (i - 1) * nlev_packed
      idx  = self%send_local_idx(i)

      is_missing = (idx == 0)
      if (.not. is_missing .and. have_fill) then
        ! The caller's own fill marker (ELM/MOSART spval) is z-invariant per
        ! column, so testing the first level is enough.
        is_missing = (fld_local(idx, 1) == lfill)
      end if

      if (is_missing) then
        do k = 1, numlev
          self%ws_send_buf(base + k) = 0.0_r8
        end do
        self%ws_send_buf(base + nlev_packed) = missing_mask
      else
        do k = 1, numlev
          self%ws_send_buf(base + k) = fld_local(idx, k)
        end do
        self%ws_send_buf(base + nlev_packed) = 1.0_r8
      end if
    end do

    call self%shared%apply_masked(self%ws_send_buf, numlev, fld_out, &
         self%mpicom, self%npes, ierr)

  end subroutine shr_horiz_remap_col_field

  !-------------------------------------------------------------------------------------------
  subroutine shr_horiz_remap_col_write(self, File, varid, fld_out, numlev, data_type)
    !
    ! Write a remapped field with PIO.  Decompositions are cached per instance
    ! and rebuilt only when the level count or data type changes.
    !
    use pio, only: file_desc_t, var_desc_t, pio_initdecomp, pio_freedecomp, &
                   pio_write_darray, PIO_OFFSET_KIND

    class(shr_horiz_remap_col_t), intent(inout) :: self
    type(file_desc_t), intent(inout) :: File
    type(var_desc_t),  intent(inout) :: varid
    real(r8),          intent(in)    :: fld_out(:,:)
    integer,           intent(in)    :: numlev
    integer,           intent(in)    :: data_type

    integer(PIO_OFFSET_KIND), allocatable :: idof(:)
    integer :: i, k, global_row, ilon, ilat, ierr
    integer :: nlat, nlon

    nlat = self%shared%nlat
    nlon = self%shared%nlon

    if (numlev <= 1) then
      if (self%iodesc_2d_valid .and. self%iodesc_2d_dtype /= data_type) then
        call pio_freedecomp(File, self%iodesc_2d)
        self%iodesc_2d_valid = .false.
      end if
      if (.not. self%iodesc_2d_valid) then
        allocate(idof(max(1, self%shared%n_b_local)))
        idof(:) = 0
        do i = 1, self%shared%n_b_local
          global_row = self%shared%row_start + i - 1
          ilon = mod(global_row - 1, nlon) + 1
          ilat = (global_row - 1) / nlon + 1
          idof(i) = int(ilon + nlon * (ilat - 1), PIO_OFFSET_KIND)
        end do
        call pio_initdecomp(self%iosystem, data_type, (/nlon, nlat/), &
             idof(1:self%shared%n_b_local), self%iodesc_2d)
        deallocate(idof)
        self%iodesc_2d_valid = .true.
        self%iodesc_2d_dtype = data_type
      end if
      call pio_write_darray(File, varid, self%iodesc_2d, fld_out(:,1), ierr)
    else
      if (self%iodesc_3d_valid .and. &
          (self%iodesc_3d_nlev /= numlev .or. self%iodesc_3d_dtype /= data_type)) then
        call pio_freedecomp(File, self%iodesc_3d)
        self%iodesc_3d_valid = .false.
      end if
      if (.not. self%iodesc_3d_valid) then
        allocate(idof(max(1, self%shared%n_b_local * numlev)))
        idof(:) = 0
        do i = 1, self%shared%n_b_local
          global_row = self%shared%row_start + i - 1
          ilon = mod(global_row - 1, nlon) + 1
          ilat = (global_row - 1) / nlon + 1
          do k = 1, numlev
            idof((k-1)*self%shared%n_b_local + i) = &
                 int(ilon + nlon*(ilat-1) + nlon*nlat*(k-1), PIO_OFFSET_KIND)
          end do
        end do
        call pio_initdecomp(self%iosystem, data_type, (/nlon, nlat, numlev/), &
             idof(1:self%shared%n_b_local*numlev), self%iodesc_3d)
        deallocate(idof)
        self%iodesc_3d_valid = .true.
        self%iodesc_3d_nlev  = numlev
        self%iodesc_3d_dtype = data_type
      end if
      call pio_write_darray(File, varid, self%iodesc_3d, fld_out, ierr)
    end if

  end subroutine shr_horiz_remap_col_write

end module shr_horiz_remap_col_mod
