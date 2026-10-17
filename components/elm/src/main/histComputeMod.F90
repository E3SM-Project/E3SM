module histComputeMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Compute side of the ELM history system: accumulates each time step's model
  ! fields into the history buffers, and normalizes / zeroes the buffers at the
  ! end of a history interval. histFileMod is the API: it owns the master field
  ! list, the tapes, the namelist, I/O and restart, and describes the active
  ! fields to this module (hist_compute_build).
  !
  ! Nothing here is ELM-specific: the subgrid structure, weights, active flags
  ! and model arrays are passed in, so the module is unit tested with a
  ! synthetic subgrid (test/hist_compute_test).
  !
  ! Layout
  ! - The buffers of all active fields on all tapes are one contiguous pool,
  !   ordered by tape, then field; each field is stored [point, level]. The
  !   tapes' hbuf/nacs are pointer views into the pool (hist_compute_views).
  ! - Field descriptors are a structure of arrays of integer codes, so the
  !   per-step loops touch no strings and no tape structures.
  ! - The four subgrid levels (topounit, landunit, column, patch) share one
  !   index space: sg = point + sg_shift(level).
  ! - A gridcell's points need not be contiguous (ELM orders landunits by type
  !   across gridcells). For gridcell averages, the points of each level are
  !   listed by segment, ascending within a segment (k_off/k_pt), with their
  !   weights, flags and scale factors stored in the same order. A segment is
  !   a tile of gridcells on the CPU: an average streams the tile's points in
  !   memory order into the tile's gridcells. With OpenACC it is one gridcell,
  !   so each vector lane averages one gridcell. Either way a gridcell adds its
  !   points in ascending order, as subgridAveMod does.
  !
  ! Parallelism
  ! - Work is cut into tiles of at most tile_len contiguous points of one level
  !   of one field. Tiles are independent: each phase is one parallel loop over
  !   tiles (OpenMP threads or OpenACC gangs) around a unit-stride sweep.
  ! - Phases: (1) snow fields rearrange their levels into a scratch pool;
  !   (2) fields averaged to gridcells reduce each gridcell's subgrid points
  !   into a gridcell scratch pool; (3) every field accumulates its values
  !   (from the model array or a scratch pool) into its buffer.
  ! - Each value is computed by one thread, in the order subgridAveMod uses (a
  !   gridcell's points ascending), so results don't depend on the number of
  !   threads and match the serial code bit for bit.
  !
  ! Accumulation rule (the same for every field): a sample is valid where the
  ! point is active (gridcell averages: some point contributed) and the value
  ! isn't spval.
  !   'I'      the sample, or spval if invalid        nacs = 1
  !   'A'      sum of valid samples (normalized: /nacs)
  !   'X','M'  max / min of valid samples
  ! For 'A', 'X' and 'M', nacs counts the valid samples; a buffer that saw no
  ! valid sample in the interval holds spval.
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  ! subgrid levels of a model array
  integer, parameter, public :: lev_grid = 1, lev_topo = 2, lev_lun = 3, lev_col = 4, lev_patch = 5
  ! time averaging
  integer, parameter, public :: avg_inst = 1, avg_mean = 2, avg_max = 3, avg_min = 4
  ! multi-layer snow fields: value used where a snow layer is absent
  integer, parameter, public :: no_snow_unset  = 0 ! not a snow field
  integer, parameter, public :: no_snow_normal = 1 ! spval (the layer doesn't count)
  integer, parameter, public :: no_snow_zero   = 2 ! 0

  ! model array registered for a history field (histFileMod: elmptr_rs/elmptr_ra).
  ! Contiguous, so the kernels get them without copies.
  type, public :: hist_src1d
     real(r8), pointer, contiguous :: ptr(:) => null()
  end type hist_src1d
  type, public :: hist_src2d
     real(r8), pointer, contiguous :: ptr(:,:) => null()
  end type hist_src2d

  ! one active field of a tape
  type, public :: hist_field_desc
     integer :: level   = lev_grid        ! subgrid level of the model array
     integer :: rank    = 1               ! 1: hist_src1d, 2: hist_src2d
     integer :: hpindex = 0               ! index of the model array
     integer :: num2d   = 1               ! number of levels
     integer :: avg     = avg_mean
     logical :: to_gridcell = .false.     ! buffer holds the gridcell average
     integer :: c2l = 1, l2g = 1          ! scale table columns (to_gridcell only)
     integer :: no_snow = no_snow_unset
  end type hist_field_desc

  ! static subgrid structure (proc bounds); built once by the API
  type, public :: hist_subgrid_desc
     integer :: begg = 1, endg = 0
     integer :: beg(lev_topo:lev_patch) = 1, end(lev_topo:lev_patch) = 0
     integer , allocatable :: grid_t(:), grid_l(:), grid_c(:), grid_p(:) ! gridcell of each point
     integer , allocatable :: lun_c(:), lun_p(:)                         ! landunit of each column / patch
     integer , allocatable :: col_p(:)                                   ! column of each patch
     real(r8), allocatable :: scale_c2l(:,:)  ! [column, c2l scale type]
     real(r8), allocatable :: scale_l2g(:,:)  ! [landunit, l2g scale type]
     logical :: extrasnowlayers = .false.     ! snow levels ordered top layer first
  end type hist_subgrid_desc

  ! subgrid state that changes during the run, read at each update
  type, public :: hist_subgrid_state
     real(r8), pointer :: wt_t(:) => null(), wt_l(:) => null(), wt_c(:) => null(), wt_p(:) => null() ! weight on gridcell
     logical , pointer :: act_t(:) => null(), act_l(:) => null(), act_c(:) => null(), act_p(:) => null()
     integer , pointer :: snl(:) => null()     ! minus the number of snow layers [begc:endc]
  end type hist_subgrid_state

  public :: hist_compute_set_subgrid ! set the subgrid structure
  public :: hist_compute_build       ! lay out buffers and work for the active fields
  public :: hist_compute_views       ! pointer views of one field's buffers
  public :: hist_compute_update      ! accumulate the current time step
  public :: hist_compute_normalize   ! divide time averages by their sample counts (one tape)
  public :: hist_compute_zero        ! reset buffers and counters (one tape)
  public :: hist_compute_update_host   ! copy one tape's buffers device -> host (OpenACC)
  public :: hist_compute_update_device ! copy all buffers host -> device (OpenACC)

  integer, parameter :: src_ptr1 = 1, src_ptr2 = 2, src_snow = 3, src_gcell = 4
  real(r8), parameter :: big = 1.e50_r8

  real(r8) :: spval = 1.e36_r8
  integer  :: tile_len = 256

  ! subgrid
  integer :: begg = 1, endg = 0, ng = 0, nsg = 0
  integer :: sg_shift(lev_grid:lev_patch) = 0
  integer , allocatable :: sg_grid(:)                    ! gridcell of each point
  real(r8), allocatable :: sg_c2l(:,:), sg_l2g(:,:)      ! [sg, scale type]: scale factors of each point
  ! points listed by segment for gridcell averages (k = position in the list)
  integer :: seg_len = 1                                 ! gridcells per segment
  integer , allocatable :: k_off(:,:)                    ! [segment, level]: start of its points
  integer , allocatable :: k_pt(:), k_g(:)               ! [k] point (level numbering), gridcell
  real(r8), allocatable :: k_c2l(:,:), k_l2g(:,:)        ! [k, scale type] scale factors
  real(r8), allocatable :: k_wt(:)                       ! [k] weight on gridcell; 0 if the point doesn't count
  integer , allocatable :: sg_col(:)                     ! column of each column / patch point
  logical , allocatable :: sg_act(:)                     ! point is active
  integer , allocatable :: sg_nsnow(:)                   ! snow layers at each column / patch point
  logical , allocatable :: all_true(:)                   ! [tile_len] mask for points not checked
  logical :: extrasnow = .false.
  integer :: lev_beg(lev_grid:lev_patch) = 1, lev_end(lev_grid:lev_patch) = 0

  ! fields
  integer :: nent = 0, ntape = 0
  integer , allocatable :: e_level(:), e_rank(:), e_hp(:), e_num2d(:), e_avg(:)
  integer , allocatable :: e_c2l(:), e_l2g(:), e_snow(:), e_src(:), e_gsrc(:)
  integer , allocatable :: e_beg(:), e_len(:), e_off(:)  ! buffer: first point, points, pool offset
  integer , allocatable :: e_sbeg(:), e_slen(:)          ! model array: first point, points
  integer , allocatable :: e_lb1(:)                      ! model array: lower bound of dimension 1
  integer , allocatable :: e_lb2(:)                      ! snow model array: lower bound of its levels
  integer , allocatable :: e_soff(:), e_goff(:)          ! scratch pool offsets
  integer , allocatable :: tape_ent(:)                   ! first field of each tape (+1 sentinel)
  real(r8), allocatable, target :: hbuf_pool(:)
  integer , allocatable, target :: nacs_pool(:)
  real(r8), allocatable :: snow_pool(:), gcell_pool(:)

  ! tiles: field, level, first and last point
  integer, allocatable :: t_ent(:), t_j(:), t_k0(:), t_k1(:)
  integer :: snow_tiles(2) = [1, 0], gcell_tiles(2) = [1, 0], acc_tiles(2) = [1, 0]
  integer, allocatable :: tape_tiles(:)                  ! first accumulate tile of each tape (+1)

contains

  !-----------------------------------------------------------------------
  subroutine hist_compute_set_subgrid(sg, ierr, errmsg)
    !
    ! The points of each gridcell at each level, and the scale factors of every
    ! point for each scale type.
    !
    type(hist_subgrid_desc), intent(in)  :: sg
    integer                , intent(out) :: ierr     ! nonzero: point with a gridcell out of range
    character(len=*)       , intent(out) :: errmsg
    integer :: lev, i, n

    begg = sg%begg; endg = sg%endg; ng = endg - begg + 1
    lev_beg(lev_grid) = begg; lev_end(lev_grid) = endg
    lev_beg(lev_topo:lev_patch) = sg%beg
    lev_end(lev_topo:lev_patch) = sg%end
    extrasnow = sg%extrasnowlayers

    n = 0
    do lev = lev_topo, lev_patch
       sg_shift(lev) = n - lev_beg(lev) + 1
       n = n + lev_end(lev) - lev_beg(lev) + 1
    end do
    nsg = n

    if (allocated(sg_grid)) deallocate(sg_grid, sg_c2l, sg_l2g, sg_col, sg_act, sg_nsnow)
    allocate(sg_grid(nsg), sg_c2l(nsg, size(sg%scale_c2l, 2)), sg_l2g(nsg, size(sg%scale_l2g, 2)))
    allocate(sg_col(nsg), sg_act(nsg), sg_nsnow(nsg))

    ierr = 0
    errmsg = ' '
    call set_grid(lev_topo,  sg%grid_t)
    call set_grid(lev_lun,   sg%grid_l)
    call set_grid(lev_col,   sg%grid_c)
    call set_grid(lev_patch, sg%grid_p)
    if (ierr /= 0) return

    ! scale factors of each point (1 where a level has none); they multiply
    ! in the order subgridAveMod uses: value*c2l*l2g*weight
    sg_c2l(:,:) = 1._r8
    sg_l2g(:,:) = 1._r8
    sg_col(:) = 0
    do i = lev_beg(lev_lun), lev_end(lev_lun)
       sg_l2g(i + sg_shift(lev_lun), :) = sg%scale_l2g(i, :)
    end do
    do i = lev_beg(lev_col), lev_end(lev_col)
       sg_c2l(i + sg_shift(lev_col), :) = sg%scale_c2l(i, :)
       sg_l2g(i + sg_shift(lev_col), :) = sg%scale_l2g(sg%lun_c(i), :)
       sg_col(i + sg_shift(lev_col)) = i
    end do
    do i = lev_beg(lev_patch), lev_end(lev_patch)
       sg_c2l(i + sg_shift(lev_patch), :) = sg%scale_c2l(sg%col_p(i), :)
       sg_l2g(i + sg_shift(lev_patch), :) = sg%scale_l2g(sg%lun_p(i), :)
       sg_col(i + sg_shift(lev_patch)) = sg%col_p(i)
    end do

  contains

    subroutine set_grid(lev, grid)
      integer, intent(in) :: lev
      integer, intent(in) :: grid(:)  ! gridcell of each point of the level (lev_beg(lev):)
      integer :: k

      do k = 1, size(grid)
         if (grid(k) < begg .or. grid(k) > endg) then
            ierr = k + lev_beg(lev) - 1
            errmsg = 'subgrid point outside the gridcell range'
            return
         end if
         sg_grid(k + lev_beg(lev) - 1 + sg_shift(lev)) = grid(k)
      end do
    end subroutine set_grid

  end subroutine hist_compute_set_subgrid

  !-----------------------------------------------------------------------
  subroutine hist_compute_build(fields, tape_first, src1, src2, spval_in, ierr, errmsg, tile)
    !
    ! Lay out the buffers (zeroed) and the work for the active fields:
    ! fields(tape_first(t) : tape_first(t+1)-1) are the fields of tape t.
    ! Call hist_compute_set_subgrid first. Existing buffers are discarded.
    !
    type(hist_field_desc), intent(in)  :: fields(:)
    integer              , intent(in)  :: tape_first(:)
    type(hist_src1d)     , intent(in)  :: src1(:)
    type(hist_src2d)     , intent(in)  :: src2(:)
    real(r8)             , intent(in)  :: spval_in
    integer              , intent(out) :: ierr
    character(len=*)     , intent(out) :: errmsg
    integer, optional    , intent(in)  :: tile     ! points per tile (tests)
    integer :: e, off, soff, goff, extent

    ierr = 0
    errmsg = ' '
    spval = spval_in
    if (present(tile)) tile_len = tile
    call release()

    nent  = size(fields)
    ntape = size(tape_first) - 1
    allocate(e_level(nent), e_rank(nent), e_hp(nent), e_num2d(nent), e_avg(nent), &
         e_c2l(nent), e_l2g(nent), e_snow(nent), e_src(nent), e_gsrc(nent), &
         e_beg(nent), e_len(nent), e_off(nent), e_sbeg(nent), e_slen(nent), e_lb1(nent), e_lb2(nent), &
         e_soff(nent), e_goff(nent), tape_ent(ntape+1))
    tape_ent(:) = tape_first(:)

    off = 0; soff = 0; goff = 0
    do e = 1, nent
       associate(fd => fields(e))
       e_level(e) = fd%level
       e_rank(e)  = fd%rank
       e_hp(e)    = fd%hpindex
       e_num2d(e) = fd%num2d
       e_avg(e)   = fd%avg
       e_c2l(e)   = fd%c2l
       e_l2g(e)   = fd%l2g
       e_snow(e)  = fd%no_snow
       e_lb2(e)   = 1
       e_sbeg(e)  = lev_beg(fd%level)
       e_slen(e)  = lev_end(fd%level) - lev_beg(fd%level) + 1

       ! the model array: its first dimension covers the level's points
       ! (from its own lower bound, which may have been reset to 1)
       if (fd%rank == 1) then
          if (.not. associated(src1(fd%hpindex)%ptr)) then
             ierr = e; errmsg = 'model array not associated'; return
          end if
          e_lb1(e) = lbound(src1(fd%hpindex)%ptr, 1)
          extent   = size(src1(fd%hpindex)%ptr, 1)
       else
          if (.not. associated(src2(fd%hpindex)%ptr)) then
             ierr = e; errmsg = 'model array not associated'; return
          end if
          e_lb1(e) = lbound(src2(fd%hpindex)%ptr, 1)
          extent   = size(src2(fd%hpindex)%ptr, 1)
          ! levels: a snow field's are counted from the array's lower bound
          ! (e.g. -nlevsno+1:0); other fields use levels 1:num2d of the array
          if (fd%no_snow /= no_snow_unset) then
             e_lb2(e) = lbound(src2(fd%hpindex)%ptr, 2)
             if (size(src2(fd%hpindex)%ptr, 2) < fd%num2d) then
                ierr = e; errmsg = 'model array has fewer levels than the field'; return
             end if
          else if (lbound(src2(fd%hpindex)%ptr, 2) > 1 .or. ubound(src2(fd%hpindex)%ptr, 2) < fd%num2d) then
             ierr = e; errmsg = 'model array does not have levels 1:num2d'; return
          end if
       end if
       if (extent < e_slen(e)) then
          ierr = e; errmsg = 'model array smaller than its subgrid level'; return
       end if
       if (fd%no_snow /= no_snow_unset .and. &
            (fd%rank /= 2 .or. (fd%level /= lev_col .and. fd%level /= lev_patch))) then
          ierr = e; errmsg = 'snow fields must be column or patch arrays with levels'; return
       end if

       ! where values come from, and where the buffer lives
       e_soff(e) = -1
       e_goff(e) = -1
       if (fd%no_snow /= no_snow_unset) then
          e_soff(e) = soff
          soff = soff + e_slen(e) * fd%num2d
          e_gsrc(e) = src_snow
       else
          e_gsrc(e) = fd%rank          ! src_ptr1 / src_ptr2
       end if
       if (fd%to_gridcell .and. fd%level /= lev_grid) then
          e_goff(e) = goff
          goff = goff + ng * fd%num2d
          e_src(e) = src_gcell
          e_beg(e) = begg
          e_len(e) = ng
       else
          e_src(e) = e_gsrc(e)
          e_beg(e) = e_sbeg(e)
          e_len(e) = e_slen(e)
       end if
       e_off(e) = off
       off = off + e_len(e) * fd%num2d
       end associate
    end do

    allocate(hbuf_pool(off), nacs_pool(off), snow_pool(soff), gcell_pool(goff), all_true(tile_len))
    all_true(:) = .true.
    hbuf_pool(:) = 0._r8
    nacs_pool(:) = 0
    call build_tiles()
    call build_segments()

    !$acc enter data copyin(hbuf_pool, nacs_pool, snow_pool, gcell_pool, &
    !$acc&   e_level, e_rank, e_hp, e_num2d, e_avg, e_c2l, e_l2g, e_snow, e_src, e_gsrc, &
    !$acc&   e_beg, e_len, e_off, e_sbeg, e_slen, e_lb1, e_lb2, e_soff, e_goff, &
    !$acc&   t_ent, t_j, t_k0, t_k1, k_off, k_pt, k_g, k_c2l, k_l2g, k_wt, sg_col, &
    !$acc&   sg_act, sg_nsnow, all_true)

  end subroutine hist_compute_build

  !-----------------------------------------------------------------------
  subroutine build_tiles()
    ! Tiles of each phase, in pool order (a tape's tiles are contiguous).
    integer :: e, t, nt, pass

    do pass = 1, 2
       nt = 0
       snow_tiles(1) = nt + 1
       do e = 1, nent
          if (e_soff(e) >= 0) call add(e, e_sbeg(e), e_slen(e))
       end do
       snow_tiles(2) = nt
       gcell_tiles(1) = nt + 1
       do e = 1, nent
          if (e_goff(e) >= 0) call add(e, begg, ng)
       end do
       gcell_tiles(2) = nt
       acc_tiles(1) = nt + 1
       if (pass == 2) tape_tiles(1) = nt + 1
       do t = 1, ntape
          do e = tape_ent(t), tape_ent(t+1) - 1
             call add(e, e_beg(e), e_len(e))
          end do
          if (pass == 2) tape_tiles(t+1) = nt + 1
       end do
       acc_tiles(2) = nt
       if (pass == 1) allocate(t_ent(nt), t_j(nt), t_k0(nt), t_k1(nt), tape_tiles(ntape+1))
    end do

  contains

    subroutine add(e, beg, len)
      integer, intent(in) :: e, beg, len
      integer :: j, k0
      do j = 1, e_num2d(e)
         do k0 = beg, beg + len - 1, tile_len
            nt = nt + 1
            if (pass == 2) then
               t_ent(nt) = e
               t_j(nt)   = j
               t_k0(nt)  = k0
               t_k1(nt)  = min(k0 + tile_len - 1, beg + len - 1)
            end if
         end do
      end do
    end subroutine add

  end subroutine build_tiles

  !-----------------------------------------------------------------------
  subroutine build_segments()
    ! Points of each level listed by segment (seg_len gridcells), ascending
    ! within a segment: a stable counting sort of the points by segment.
    integer :: lev, nseg, i, k, sgi, seg
    integer, allocatable :: next(:)

#ifdef _OPENACC
    seg_len = 1          ! one gridcell per vector lane
#else
    seg_len = tile_len   ! a tile's points are streamed in memory order
#endif
    nseg = (ng + seg_len - 1) / seg_len
    if (allocated(k_off)) deallocate(k_off, k_pt, k_g, k_c2l, k_l2g, k_wt)
    allocate(k_off(nseg+1, lev_topo:lev_patch), k_pt(nsg), k_g(nsg), &
         k_c2l(nsg, size(sg_c2l, 2)), k_l2g(nsg, size(sg_l2g, 2)), k_wt(nsg))
    allocate(next(nseg+1))
    k = 0
    do lev = lev_topo, lev_patch
       next(:) = 0
       do i = lev_beg(lev), lev_end(lev)
          seg = (sg_grid(i + sg_shift(lev)) - begg) / seg_len + 1
          next(seg) = next(seg) + 1
       end do
       k_off(1, lev) = k + 1
       do seg = 1, nseg
          k_off(seg+1, lev) = k_off(seg, lev) + next(seg)
          next(seg) = k_off(seg, lev)
       end do
       do i = lev_beg(lev), lev_end(lev)
          sgi = i + sg_shift(lev)
          seg = (sg_grid(sgi) - begg) / seg_len + 1
          k_pt(next(seg))  = i
          k_g(next(seg))   = sg_grid(sgi)
          k_c2l(next(seg), :) = sg_c2l(sgi, :)
          k_l2g(next(seg), :) = sg_l2g(sgi, :)
          next(seg) = next(seg) + 1
       end do
       k = k_off(nseg+1, lev) - 1
    end do
  end subroutine build_segments

  !-----------------------------------------------------------------------
  subroutine release()
    if (.not. allocated(e_level)) return
    !$acc exit data delete(hbuf_pool, nacs_pool, snow_pool, gcell_pool, &
    !$acc&   e_level, e_rank, e_hp, e_num2d, e_avg, e_c2l, e_l2g, e_snow, e_src, e_gsrc, &
    !$acc&   e_beg, e_len, e_off, e_sbeg, e_slen, e_lb1, e_lb2, e_soff, e_goff, &
    !$acc&   t_ent, t_j, t_k0, t_k1, k_off, k_pt, k_g, k_c2l, k_l2g, k_wt, sg_col, &
    !$acc&   sg_act, sg_nsnow, all_true)
    deallocate(e_level, e_rank, e_hp, e_num2d, e_avg, e_c2l, e_l2g, e_snow, e_src, e_gsrc, &
         e_beg, e_len, e_off, e_sbeg, e_slen, e_lb1, e_lb2, e_soff, e_goff, tape_ent, &
         hbuf_pool, nacs_pool, snow_pool, gcell_pool, t_ent, t_j, t_k0, t_k1, tape_tiles, all_true)
  end subroutine release

  !-----------------------------------------------------------------------
  subroutine hist_compute_views(e, hbuf, nacs)
    ! Buffer and counter of field e (pool order) as [point, level] arrays.
    integer          , intent(in) :: e
    real(r8), pointer, intent(out) :: hbuf(:,:)
    integer , pointer, intent(out) :: nacs(:,:)
    integer :: n

    n = e_len(e) * e_num2d(e)
    hbuf(e_beg(e):e_beg(e)+e_len(e)-1, 1:e_num2d(e)) => hbuf_pool(e_off(e)+1 : e_off(e)+n)
    nacs(e_beg(e):e_beg(e)+e_len(e)-1, 1:e_num2d(e)) => nacs_pool(e_off(e)+1 : e_off(e)+n)
  end subroutine hist_compute_views

  !-----------------------------------------------------------------------
  subroutine hist_compute_update(src1, src2, state, bad_gridcell)
    !
    ! Accumulate the current values of all active fields into their buffers.
    ! bad_gridcell is a gridcell whose weights sum to more than 1 (0 if none).
    !
    type(hist_src1d)        , intent(in)  :: src1(:)
    type(hist_src2d)        , intent(in)  :: src2(:)
    type(hist_subgrid_state), intent(in)  :: state
    integer                 , intent(out) :: bad_gridcell
    integer :: it, e, j, k, n, hp, lb, sb, se, nl, nn, o, so, sh, lev, k0, k1, bad
    real(r8) :: v

    call gather_state(state)

    ! (1) snow fields: absent layers take the no-snow value. Layers are
    ! ordered bottom-up (absent layers first) unless extrasnow, where the top
    ! layer comes first (absent layers last).
#ifdef _OPENACC
    !$acc parallel loop gang default(present) private(e, j, hp, lb, sb, nl, so, sh, v)
#else
    !$omp parallel do schedule(dynamic, 4) private(e, j, k, hp, lb, sb, nl, nn, so, sh, v)
#endif
    do it = snow_tiles(1), snow_tiles(2)
       e = t_ent(it); j = t_j(it)
       hp = e_hp(e); lb = e_lb1(e); sb = e_sbeg(e); nl = e_num2d(e)
       so = e_soff(e) + (j-1)*e_slen(e) - sb + 1
       sh = sg_shift(e_level(e))
       v = spval
       if (e_snow(e) == no_snow_zero) v = 0._r8
       !$acc loop vector private(nn)
       do k = t_k0(it), t_k1(it)
          nn = nl - sg_nsnow(k + sh)           ! absent layers
          if (.not. extrasnow) then
             if (j <= nn) then
                snow_pool(so + k) = v
             else
                snow_pool(so + k) = src2(hp)%ptr(k - sb + lb, j - 1 + e_lb2(e))
             end if
          else
             if (j > nl - nn) then
                snow_pool(so + k) = v
             else
                snow_pool(so + k) = src2(hp)%ptr(k - sb + lb, j + nn - 1 + e_lb2(e))
             end if
          end if
       end do
    end do

    ! (2) gridcell averages: weighted mean over the gridcell's valid points
    bad_gridcell = 0
#ifdef _OPENACC
    !$acc parallel loop gang default(present) reduction(max:bad_gridcell) &
    !$acc&   private(e, j, lev, hp, lb, sb, se, so, o, bad)
#else
    !$omp parallel do schedule(dynamic, 4) reduction(max:bad_gridcell) &
    !$omp&   private(e, j, lev, hp, lb, sb, se, so, o, bad)
#endif
    do it = gcell_tiles(1), gcell_tiles(2)
       e = t_ent(it); j = t_j(it); lev = e_level(e)
       hp = e_hp(e); lb = e_lb1(e); sb = e_sbeg(e); se = sb + e_slen(e) - 1
       o  = e_goff(e) + (j-1)*ng - begg + 1
       select case (e_gsrc(e))
       case (src_ptr1)
          call average_tile(t_k0(it), t_k1(it), lev, sb, src1(hp)%ptr(lb:lb+se-sb), &
               e_c2l(e), e_l2g(e), gcell_pool(o + t_k0(it) : o + t_k1(it)), bad)
       case (src_ptr2)
          call average_tile(t_k0(it), t_k1(it), lev, sb, src2(hp)%ptr(lb:lb+se-sb, j), &
               e_c2l(e), e_l2g(e), gcell_pool(o + t_k0(it) : o + t_k1(it)), bad)
       case default
          so = e_soff(e) + (j-1)*e_slen(e)
          call average_tile(t_k0(it), t_k1(it), lev, sb, snow_pool(so+1 : so+e_slen(e)), &
               e_c2l(e), e_l2g(e), gcell_pool(o + t_k0(it) : o + t_k1(it)), bad)
       end select
       bad_gridcell = max(bad_gridcell, bad)
    end do

    ! (3) accumulate every field into its buffer
#ifdef _OPENACC
    !$acc parallel loop gang default(present) private(e, j, hp, lb, sb, sh, so, o, k0, k1, n)
#else
    !$omp parallel do schedule(dynamic, 4) private(e, j, hp, lb, sb, sh, so, o, k0, k1, n)
#endif
    do it = acc_tiles(1), acc_tiles(2)
       e = t_ent(it); j = t_j(it); k0 = t_k0(it); k1 = t_k1(it); n = k1 - k0 + 1
       hp = e_hp(e); lb = e_lb1(e); sb = e_beg(e)
       sh = sg_shift(e_level(e))
       o  = e_off(e) + (j-1)*e_len(e) - sb + 1   ! hbuf_pool(o + k) is point k
       select case (e_src(e))
       case (src_gcell)
          so = e_goff(e) + (j-1)*ng - begg + 1
          call accumulate_tile(e_avg(e), n, gcell_pool(so+k0 : so+k1), all_true(1:n), &
               spval, hbuf_pool(o+k0 : o+k1), nacs_pool(o+k0 : o+k1))
       case (src_snow)
          so = e_soff(e) + (j-1)*e_slen(e) - sb + 1
          call accumulate_tile(e_avg(e), n, snow_pool(so+k0 : so+k1), sg_act(k0+sh : k1+sh), &
               spval, hbuf_pool(o+k0 : o+k1), nacs_pool(o+k0 : o+k1))
       case (src_ptr1)
          if (e_level(e) == lev_grid) then
             call accumulate_tile(e_avg(e), n, src1(hp)%ptr(k0-sb+lb : k1-sb+lb), all_true(1:n), &
                  spval, hbuf_pool(o+k0 : o+k1), nacs_pool(o+k0 : o+k1))
          else
             call accumulate_tile(e_avg(e), n, src1(hp)%ptr(k0-sb+lb : k1-sb+lb), sg_act(k0+sh : k1+sh), &
                  spval, hbuf_pool(o+k0 : o+k1), nacs_pool(o+k0 : o+k1))
          end if
       case (src_ptr2)
          if (e_level(e) == lev_grid) then
             call accumulate_tile(e_avg(e), n, src2(hp)%ptr(k0-sb+lb : k1-sb+lb, j), all_true(1:n), &
                  spval, hbuf_pool(o+k0 : o+k1), nacs_pool(o+k0 : o+k1))
          else
             call accumulate_tile(e_avg(e), n, src2(hp)%ptr(k0-sb+lb : k1-sb+lb, j), sg_act(k0+sh : k1+sh), &
                  spval, hbuf_pool(o+k0 : o+k1), nacs_pool(o+k0 : o+k1))
          end if
       end select
    end do

  end subroutine hist_compute_update

  !-----------------------------------------------------------------------
  subroutine average_tile(g0, g1, lev, sb, v, ic2l, il2g, garr, bad)
    !
    ! Gridcell averages of one level of a field, for gridcells g0:g1:
    ! garr(g) = sum(v*c2l*l2g*wt) / sum(wt) over g's points that count, in
    ! ascending point order; spval if none does. bad: a gridcell whose
    ! weights sum to more than 1 (else 0). Scale type 1 is 'unity': a factor
    ! of exactly 1 is skipped (same result, half the memory traffic). Points
    ! that don't count are selected away rather than branched around, so the
    ! cost doesn't depend on the pattern of missing values (the products of
    ! spval stay far from overflow).
    !
    !$acc routine vector
    integer , intent(in)  :: g0, g1
    integer , intent(in)  :: lev           ! subgrid level of the values
    integer , intent(in)  :: sb            ! first point of the level
    real(r8), intent(in), contiguous :: v(sb:)  ! the field's values at the level's points (sb:se)
    integer , intent(in)  :: ic2l, il2g    ! scale types
    real(r8), intent(out) :: garr(g0:g1)
    integer , intent(out) :: bad
    integer  :: g, k, seg
    logical  :: use_c2l, use_l2g, ok
    real(r8) :: x, f
#ifdef _OPENACC
    real(r8) :: s, sw
#else
    real(r8) :: s(g0:g1), sw(g0:g1)
#endif

    use_c2l = ic2l /= 1    ! factors that aren't unity
    use_l2g = il2g /= 1
    bad = 0
#ifdef _OPENACC
    ! one gridcell per vector lane (segments are single gridcells)
    !$acc loop vector private(s, sw, x, f, ok, seg) reduction(max:bad)
    do g = g0, g1
       seg = g - begg + 1
       s  = spval
       sw = 0._r8
       do k = k_off(seg, lev), k_off(seg+1, lev) - 1
          x  = v(k_pt(k))
          ok = k_wt(k) /= 0._r8 .and. x /= spval
          f  = x
          if (use_c2l) then
             ok = ok .and. k_c2l(k, ic2l) /= spval
             f  = f * k_c2l(k, ic2l)
          end if
          if (use_l2g) then
             ok = ok .and. k_l2g(k, il2g) /= spval
             f  = f * k_l2g(k, il2g)
          end if
          s  = merge(merge(0._r8, s, sw == 0._r8) + f * k_wt(k), s, ok)
          sw = merge(sw + k_wt(k), sw, ok)
       end do
       if (sw > 1.0_r8 + 1.e-6_r8) then
          bad = max(bad, g)
       else if (sw /= 0._r8) then
          s = s / sw
       end if
       garr(g) = s
    end do
#else
    ! stream the tile's points (one segment) in memory order
    seg = (g0 - begg) / seg_len + 1
    s(:)  = spval
    sw(:) = 0._r8
    do k = k_off(seg, lev), k_off(seg+1, lev) - 1
       x  = v(k_pt(k))
       ok = k_wt(k) /= 0._r8 .and. x /= spval
       f  = x
       if (use_c2l) then
          ok = ok .and. k_c2l(k, ic2l) /= spval
          f  = f * k_c2l(k, ic2l)
       end if
       if (use_l2g) then
          ok = ok .and. k_l2g(k, il2g) /= spval
          f  = f * k_l2g(k, il2g)
       end if
       g = k_g(k)
       s(g)  = merge(merge(0._r8, s(g), sw(g) == 0._r8) + f * k_wt(k), s(g), ok)
       sw(g) = merge(sw(g) + k_wt(k), sw(g), ok)
    end do
    do g = g0, g1
       if (sw(g) > 1.0_r8 + 1.e-6_r8) then
          bad = max(bad, g)
       else if (sw(g) /= 0._r8) then
          s(g) = s(g) / sw(g)
       end if
       garr(g) = s(g)
    end do
#endif


  end subroutine average_tile

  !-----------------------------------------------------------------------
  subroutine accumulate_tile(avg, n, v, act, missing, hbuf, nacs)
    !
    ! Accumulate n consecutive samples into their buffers (see the rule in the
    ! module description). A sample is valid where act is true and v /= missing.
    ! Selects rather than branches, so the loops vectorize.
    !
    !$acc routine vector
    integer , intent(in)    :: avg, n
    real(r8), intent(in)   , contiguous :: v(:)     ! (n)
    logical , intent(in)   , contiguous :: act(:)
    real(r8), intent(in)    :: missing
    real(r8), intent(inout), contiguous :: hbuf(:)
    integer , intent(inout), contiguous :: nacs(:)
    integer :: k
    logical :: valid

    select case (avg)
    case (avg_inst)
       !$acc loop vector
       do k = 1, n
          hbuf(k) = merge(v(k), missing, act(k) .and. v(k) /= missing)
          nacs(k) = 1
       end do
    case (avg_mean)
       !$acc loop vector private(valid)
       do k = 1, n
          valid = act(k) .and. v(k) /= missing
          hbuf(k) = merge(merge(0._r8, hbuf(k), nacs(k) == 0) + v(k), &
                          merge(missing, hbuf(k), nacs(k) == 0), valid)
          nacs(k) = nacs(k) + merge(1, 0, valid)
       end do
    case (avg_max)
       !$acc loop vector private(valid)
       do k = 1, n
          valid = act(k) .and. v(k) /= missing
          hbuf(k) = merge(max(merge(-big, hbuf(k), nacs(k) == 0), v(k)), &
                          merge(missing, hbuf(k), nacs(k) == 0), valid)
          nacs(k) = nacs(k) + merge(1, 0, valid)
       end do
    case (avg_min)
       !$acc loop vector private(valid)
       do k = 1, n
          valid = act(k) .and. v(k) /= missing
          hbuf(k) = merge(min(merge(big, hbuf(k), nacs(k) == 0), v(k)), &
                          merge(missing, hbuf(k), nacs(k) == 0), valid)
          nacs(k) = nacs(k) + merge(1, 0, valid)
       end do
    end select
  end subroutine accumulate_tile


  !-----------------------------------------------------------------------
  subroutine gather_state(state)
    ! Active flags and snow layers of all points (sg order), and the weights
    ! of the points and whether they count in gridcell averages (segment order).
    type(hist_subgrid_state), intent(in) :: state

    call gather(lev_topo,  state%wt_t, state%act_t, .false.)
    call gather(lev_lun,   state%wt_l, state%act_l, .true.)
    call gather(lev_col,   state%wt_c, state%act_c, .true.)
    call gather(lev_patch, state%wt_p, state%act_p, .true.)
    !$acc update device(k_wt, sg_act, sg_nsnow)

  contains

    subroutine gather(lev, wt, act, check_active)
      integer , intent(in) :: lev
      real(r8), intent(in) :: wt(lev_beg(lev):)
      logical , intent(in) :: act(lev_beg(lev):)
      logical , intent(in) :: check_active  ! gridcell averages skip inactive points
      integer :: i, n

      !$omp parallel do private(i)
      do n = k_off(1, lev), k_off(size(k_off, 1), lev) - 1
         i = k_pt(n)
         ! (weight /= 0 and active) is the test subgridAveMod uses
         k_wt(n) = merge(wt(i), 0._r8, act(i) .or. .not. check_active)
      end do
      !$omp parallel do private(n)
      do i = lev_beg(lev), lev_end(lev)
         n = i + sg_shift(lev)
         sg_act(n) = act(i)
         if (lev >= lev_col .and. associated(state%snl)) then
            sg_nsnow(n) = abs(state%snl(sg_col(n)))
         else
            sg_nsnow(n) = 0
         end if
      end do
    end subroutine gather

  end subroutine gather_state

  !-----------------------------------------------------------------------
  subroutine hist_compute_normalize(t)
    ! Divide the time averages of tape t by their number of samples.
    integer, intent(in) :: t
    integer :: it, e, j, k, o

#ifdef _OPENACC
    !$acc parallel loop gang default(present) private(e, j, o)
#else
    !$omp parallel do private(e, j, k, o)
#endif
    do it = tape_tiles(t), tape_tiles(t+1) - 1
       e = t_ent(it)
       if (e_avg(e) /= avg_mean) cycle
       j = t_j(it)
       o = e_off(e) + (j-1)*e_len(e) - e_beg(e) + 1
       !$acc loop vector
       do k = t_k0(it), t_k1(it)
          if (nacs_pool(o + k) /= 0) hbuf_pool(o + k) = hbuf_pool(o + k) / real(nacs_pool(o + k), r8)
       end do
    end do
  end subroutine hist_compute_normalize

  !-----------------------------------------------------------------------
  subroutine hist_compute_zero(t)
    ! Reset the buffers and counters of tape t.
    integer, intent(in) :: t
    integer :: n, n0, n1

    if (tape_ent(t+1) <= tape_ent(t)) return
    n0 = e_off(tape_ent(t)) + 1
    n1 = e_off(tape_ent(t+1) - 1) + e_len(tape_ent(t+1) - 1) * e_num2d(tape_ent(t+1) - 1)
#ifdef _OPENACC
    !$acc parallel loop default(present)
#else
    !$omp parallel do
#endif
    do n = n0, n1
       hbuf_pool(n) = 0._r8
       nacs_pool(n) = 0
    end do
  end subroutine hist_compute_zero

  !-----------------------------------------------------------------------
  subroutine hist_compute_update_host(t)
    ! Before tape t's buffers are written (OpenACC builds; no-op otherwise).
    integer, intent(in) :: t
    integer :: n0, n1

    if (tape_ent(t+1) <= tape_ent(t)) return
    n0 = e_off(tape_ent(t)) + 1
    n1 = e_off(tape_ent(t+1) - 1) + e_len(tape_ent(t+1) - 1) * e_num2d(tape_ent(t+1) - 1)
    !$acc update self(hbuf_pool(n0:n1), nacs_pool(n0:n1))
  end subroutine hist_compute_update_host

  !-----------------------------------------------------------------------
  subroutine hist_compute_update_device()
    ! After the buffers were read on the host (restart; OpenACC builds only).
    !$acc update device(hbuf_pool, nacs_pool)
  end subroutine hist_compute_update_device

end module histComputeMod
