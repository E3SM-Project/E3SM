module subgridAveViewMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! The subgrid structure that averages to gridcells need, passed around as an
  ! argument (subgrid_ave_view) instead of read from ELM's global instances
  ! (grc_pp, top_pp, lun_pp, col_pp, veg_pp), and the averages to gridcells
  ! that subgridAveMod (p2g, c2g, l2g, t2g) and the history (histComputeMod)
  ! share. Nothing here uses ELM's subgrid types, so it is unit tested with a
  ! synthetic subgrid (main/test/hist_compute_test).
  !
  ! Order of the sums: a gridcell's landunits are listed in ascending index
  ! order (grc_lun_ptr/grc_lun_idx, the one index built here), and the columns
  ! and patches of a landunit are contiguous (coli:colf, pfti:pftf), so an
  ! average adds a gridcell's points in ascending index order, the order of the
  ! serial loops it replaces: results are bit-for-bit with them, whatever order
  ! ELM creates landunit types in. A gridcell's topounits are contiguous
  ! (topi:topf).
  !
  ! Each average handles a range of gridcells with the gridcell loop inside
  ! (OpenACC: a vector routine), so its caller makes one call per field level
  ! instead of one per gridcell. The averages don't check the sums of weights:
  ! reweight_wrapup (check_weights) stops the run if the active weights of a
  ! gridcell don't sum to 1, every time they change. They are pure. Values are passed as assumed-size arrays whose
  ! first element is the first point of their level (e.g. parr(1) is patch
  ! view%begp), so a model array, or a contiguous section of one, is passed
  ! without a copy.
  !
  ! In ELM, subgridAveMod builds the view (elm_subgrid_view): everything but
  ! the landunit index points at the subgrid types' arrays, so changes of the
  ! weights and active flags are seen without a refresh. The view's pointers
  ! are contiguous (their targets must be), so the loops need no strides.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  implicit none
  private
  !
  ! !PUBLIC DATA:
  real(r8), parameter, public :: spval = 1.e36_r8              ! missing value (= elm_varcon's spval)
  integer , parameter, public :: scale_unity = 1               ! the 'unity' c2l/l2g scale type id (= subgridAveMod's unity)

  type, public :: subgrid_ave_view
     integer :: begg = 1, endg = 0, begt = 1, endt = 0, begl = 1, endl = 0
     integer :: begc = 1, endc = 0, begp = 1, endp = 0
     ! landunits of each gridcell, in ascending order: grc_lun_idx(grc_lun_ptr(g):grc_lun_ptr(g+1)-1)
     integer , pointer, contiguous :: grc_lun_ptr(:) => null(), grc_lun_idx(:) => null()
     ! contiguous children: topounits of a gridcell, columns and patches of a landunit
     integer , pointer, contiguous :: topi(:) => null(), topf(:) => null()
     integer , pointer, contiguous :: coli(:) => null(), colf(:) => null()
     integer , pointer, contiguous :: pfti(:) => null(), pftf(:) => null()
     ! column of each patch, type of each landunit (lun_pp%itype)
     integer , pointer, contiguous :: col_of_pft(:) => null(), lun_type(:) => null()
     ! weight of each point on its gridcell
     real(r8), pointer, contiguous :: wtgcell_t(:) => null(), wtgcell_l(:) => null()
     real(r8), pointer, contiguous :: wtgcell_c(:) => null(), wtgcell_p(:) => null()
     ! point is active
     logical , pointer, contiguous :: active_t(:) => null(), active_l(:) => null()
     logical , pointer, contiguous :: active_c(:) => null(), active_p(:) => null()
     ! scale factors: scale_c2l(column, c2l scale type) and
     ! scale_l2g_lookup(l2g scale type, landunit type); spval excludes a point
     real(r8), pointer, contiguous :: scale_c2l(:,:) => null(), scale_l2g_lookup(:,:) => null()
     ! minus the number of snow layers of each column
     integer , pointer, contiguous :: snl(:) => null()
  end type subgrid_ave_view
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: build_child_map          ! children of each parent, in ascending order
  public :: subgrid_ave_view_attach  ! put a view on the device (OpenACC)
  public :: subgrid_ave_view_detach  ! remove a view from the device (OpenACC)
  public :: p2g_range                ! patch averages of a range of gridcells
  public :: c2g_range                ! column averages of a range of gridcells
  public :: p2g_range_unity          ! p2g_range for the unity scale types (no lookups)
  public :: c2g_range_unity          ! c2g_range for the unity scale types (no lookups)
  public :: l2g_range                ! landunit averages of a range of gridcells
  public :: t2g_range                ! topounit averages of a range of gridcells
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  pure subroutine build_child_map(begpar, endpar, begch, endch, parent, ptr, idx, bad_child)
    !
    ! !DESCRIPTION:
    ! Counting sort of the children by parent, keeping ascending child order.
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: begpar, endpar, begch, endch
    integer, intent(in)  :: parent(begch:endch)  ! parent of each child
    integer, intent(out) :: ptr(begpar:endpar+1) ! children of k: idx(ptr(k):ptr(k+1)-1)
    integer, intent(out) :: idx(begch:endch)
    integer, intent(out) :: bad_child            ! a child whose parent is out of range, or 0
    !
    ! !LOCAL VARIABLES:
    integer :: i, k
    integer :: next(begpar:endpar)
    !-----------------------------------------------------------------------

    bad_child = 0
    next(:) = 0
    do i = begch, endch
       if (parent(i) < begpar .or. parent(i) > endpar) then
          bad_child = i
          return
       end if
       next(parent(i)) = next(parent(i)) + 1
    end do
    ptr(begpar) = begch
    do k = begpar, endpar
       ptr(k+1) = ptr(k) + next(k)
    end do
    next(:) = ptr(begpar:endpar)
    do i = begch, endch
       idx(next(parent(i))) = i
       next(parent(i)) = next(parent(i)) + 1
    end do

  end subroutine build_child_map

  !-----------------------------------------------------------------------
  subroutine subgrid_ave_view_attach(view)
    !
    ! !DESCRIPTION:
    ! Copy the view to the device and attach its pointers to the device copies
    ! of their targets, which must already be on the device (a pointer whose
    ! target isn't present is left pointing at host memory).
    !
    ! !ARGUMENTS:
    type(subgrid_ave_view), intent(inout) :: view
    !-----------------------------------------------------------------------

    !$acc enter data copyin(view)
    !$acc enter data attach(view%grc_lun_ptr, view%grc_lun_idx, view%topi, view%topf, &
    !$acc    view%coli, view%colf, view%pfti, view%pftf, view%col_of_pft, view%lun_type, &
    !$acc    view%wtgcell_t, view%wtgcell_l, view%wtgcell_c, view%wtgcell_p, &
    !$acc    view%active_t, view%active_l, view%active_c, view%active_p, &
    !$acc    view%scale_c2l, view%scale_l2g_lookup, view%snl)

  end subroutine subgrid_ave_view_attach

  !-----------------------------------------------------------------------
  subroutine subgrid_ave_view_detach(view)
    !
    ! !DESCRIPTION:
    ! Undo subgrid_ave_view_attach.
    !
    ! !ARGUMENTS:
    type(subgrid_ave_view), intent(inout) :: view
    !-----------------------------------------------------------------------

    !$acc exit data detach(view%grc_lun_ptr, view%grc_lun_idx, view%topi, view%topf, &
    !$acc    view%coli, view%colf, view%pfti, view%pftf, view%col_of_pft, view%lun_type, &
    !$acc    view%wtgcell_t, view%wtgcell_l, view%wtgcell_c, view%wtgcell_p, &
    !$acc    view%active_t, view%active_l, view%active_c, view%active_p, &
    !$acc    view%scale_c2l, view%scale_l2g_lookup, view%snl)
    !$acc exit data delete(view)

  end subroutine subgrid_ave_view_detach

  !-----------------------------------------------------------------------
  pure subroutine p2g_range(view, g0, g1, parr, beg, c2l_scale_type, l2g_scale_type, garr)
    !
    ! !DESCRIPTION:
    ! For gridcells g0:g1, the average of each gridcell's active patches with
    ! nonzero weight whose value and scale factors aren't spval (spval if
    ! there are none).
    !
    !$acc routine vector
    ! !ARGUMENTS:
    type(subgrid_ave_view), intent(in) :: view
    integer , intent(in)  :: g0, g1
    real(r8), intent(in)  :: parr(*)          ! patch values; parr(1) is patch beg
    integer , intent(in)  :: beg              ! index of parr(1)
    integer , intent(in)  :: c2l_scale_type
    integer , intent(in)  :: l2g_scale_type
    real(r8), intent(out) :: garr(*)          ! garr(1) is gridcell g0
    !
    ! !LOCAL VARIABLES:
    integer  :: g, k, l, p
    real(r8) :: value, sum_arr, sumwt, scale_c2l, scale_l2g
    !-----------------------------------------------------------------------

    !$acc loop vector private(sum_arr, sumwt)
    do g = g0, g1
       sum_arr = 0._r8
       sumwt = 0._r8
       do k = view%grc_lun_ptr(g), view%grc_lun_ptr(g+1)-1
          l = view%grc_lun_idx(k)
          scale_l2g = view%scale_l2g_lookup(l2g_scale_type, view%lun_type(l))
          if (scale_l2g == spval) cycle
          do p = view%pfti(l), view%pftf(l)
             if (view%active_p(p) .and. view%wtgcell_p(p) /= 0._r8) then
                scale_c2l = view%scale_c2l(view%col_of_pft(p), c2l_scale_type)
                value = parr(p - beg + 1)
                if (value /= spval .and. scale_c2l /= spval) then
                   sum_arr = sum_arr + value * scale_c2l * scale_l2g * view%wtgcell_p(p)
                   sumwt = sumwt + view%wtgcell_p(p)
                end if
             end if
          end do
       end do
       garr(g - g0 + 1) = weighted_average(sum_arr, sumwt)
    end do

  end subroutine p2g_range

  !-----------------------------------------------------------------------
  pure subroutine c2g_range(view, g0, g1, carr, beg, c2l_scale_type, l2g_scale_type, garr)
    !
    ! !DESCRIPTION:
    ! For gridcells g0:g1, the average of each gridcell's active columns with
    ! nonzero weight whose value and scale factors aren't spval (spval if
    ! there are none).
    !
    !$acc routine vector
    ! !ARGUMENTS:
    type(subgrid_ave_view), intent(in) :: view
    integer , intent(in)  :: g0, g1
    real(r8), intent(in)  :: carr(*)          ! column values; carr(1) is column beg
    integer , intent(in)  :: beg              ! index of carr(1)
    integer , intent(in)  :: c2l_scale_type
    integer , intent(in)  :: l2g_scale_type
    real(r8), intent(out) :: garr(*)          ! garr(1) is gridcell g0
    !
    ! !LOCAL VARIABLES:
    integer  :: g, k, l, c
    real(r8) :: value, sum_arr, sumwt, scale_c2l, scale_l2g
    !-----------------------------------------------------------------------

    !$acc loop vector private(sum_arr, sumwt)
    do g = g0, g1
       sum_arr = 0._r8
       sumwt = 0._r8
       do k = view%grc_lun_ptr(g), view%grc_lun_ptr(g+1)-1
          l = view%grc_lun_idx(k)
          scale_l2g = view%scale_l2g_lookup(l2g_scale_type, view%lun_type(l))
          if (scale_l2g == spval) cycle
          do c = view%coli(l), view%colf(l)
             if (view%active_c(c) .and. view%wtgcell_c(c) /= 0._r8) then
                scale_c2l = view%scale_c2l(c, c2l_scale_type)
                value = carr(c - beg + 1)
                if (value /= spval .and. scale_c2l /= spval) then
                   sum_arr = sum_arr + value * scale_c2l * scale_l2g * view%wtgcell_c(c)
                   sumwt = sumwt + view%wtgcell_c(c)
                end if
             end if
          end do
       end do
       garr(g - g0 + 1) = weighted_average(sum_arr, sumwt)
    end do

  end subroutine c2g_range

  !-----------------------------------------------------------------------
  pure subroutine p2g_range_unity(view, g0, g1, parr, beg, garr)
    !
    ! !DESCRIPTION:
    ! p2g_range for the 'unity' c2l and l2g scale types (the default of every
    ! history field): both scale factors are exactly 1, so the lookups, their
    ! spval checks and the two multiplications are dropped. The result is bit
    ! for bit that of p2g_range with scale_unity.
    !
    !$acc routine vector
    ! !ARGUMENTS:
    type(subgrid_ave_view), intent(in) :: view
    integer , intent(in)  :: g0, g1
    real(r8), intent(in)  :: parr(*)          ! patch values; parr(1) is patch beg
    integer , intent(in)  :: beg              ! index of parr(1)
    real(r8), intent(out) :: garr(*)          ! garr(1) is gridcell g0
    !
    ! !LOCAL VARIABLES:
    integer  :: g, k, l, p
    real(r8) :: value, sum_arr, sumwt
    !-----------------------------------------------------------------------

    !$acc loop vector private(sum_arr, sumwt)
    do g = g0, g1
       sum_arr = 0._r8
       sumwt = 0._r8
       do k = view%grc_lun_ptr(g), view%grc_lun_ptr(g+1)-1
          l = view%grc_lun_idx(k)
          do p = view%pfti(l), view%pftf(l)
             if (view%active_p(p) .and. view%wtgcell_p(p) /= 0._r8) then
                value = parr(p - beg + 1)
                if (value /= spval) then
                   sum_arr = sum_arr + value * view%wtgcell_p(p)
                   sumwt = sumwt + view%wtgcell_p(p)
                end if
             end if
          end do
       end do
       garr(g - g0 + 1) = weighted_average(sum_arr, sumwt)
    end do

  end subroutine p2g_range_unity

  !-----------------------------------------------------------------------
  pure subroutine c2g_range_unity(view, g0, g1, carr, beg, garr)
    !
    ! !DESCRIPTION:
    ! c2g_range for the 'unity' c2l and l2g scale types; see p2g_range_unity.
    !
    !$acc routine vector
    ! !ARGUMENTS:
    type(subgrid_ave_view), intent(in) :: view
    integer , intent(in)  :: g0, g1
    real(r8), intent(in)  :: carr(*)          ! column values; carr(1) is column beg
    integer , intent(in)  :: beg              ! index of carr(1)
    real(r8), intent(out) :: garr(*)          ! garr(1) is gridcell g0
    !
    ! !LOCAL VARIABLES:
    integer  :: g, k, l, c
    real(r8) :: value, sum_arr, sumwt
    !-----------------------------------------------------------------------

    !$acc loop vector private(sum_arr, sumwt)
    do g = g0, g1
       sum_arr = 0._r8
       sumwt = 0._r8
       do k = view%grc_lun_ptr(g), view%grc_lun_ptr(g+1)-1
          l = view%grc_lun_idx(k)
          do c = view%coli(l), view%colf(l)
             if (view%active_c(c) .and. view%wtgcell_c(c) /= 0._r8) then
                value = carr(c - beg + 1)
                if (value /= spval) then
                   sum_arr = sum_arr + value * view%wtgcell_c(c)
                   sumwt = sumwt + view%wtgcell_c(c)
                end if
             end if
          end do
       end do
       garr(g - g0 + 1) = weighted_average(sum_arr, sumwt)
    end do

  end subroutine c2g_range_unity

  !-----------------------------------------------------------------------
  pure subroutine l2g_range(view, g0, g1, larr, beg, l2g_scale_type, garr)
    !
    ! !DESCRIPTION:
    ! For gridcells g0:g1, the average of each gridcell's active landunits with
    ! nonzero weight whose value and scale factor aren't spval (spval if there
    ! are none).
    !
    !$acc routine vector
    ! !ARGUMENTS:
    type(subgrid_ave_view), intent(in) :: view
    integer , intent(in)  :: g0, g1
    real(r8), intent(in)  :: larr(*)          ! landunit values; larr(1) is landunit beg
    integer , intent(in)  :: beg              ! index of larr(1)
    integer , intent(in)  :: l2g_scale_type
    real(r8), intent(out) :: garr(*)          ! garr(1) is gridcell g0
    !
    ! !LOCAL VARIABLES:
    integer  :: g, k, l
    real(r8) :: value, sum_arr, sumwt, scale_l2g
    !-----------------------------------------------------------------------

    !$acc loop vector private(sum_arr, sumwt)
    do g = g0, g1
       sum_arr = 0._r8
       sumwt = 0._r8
       do k = view%grc_lun_ptr(g), view%grc_lun_ptr(g+1)-1
          l = view%grc_lun_idx(k)
          if (view%active_l(l) .and. view%wtgcell_l(l) /= 0._r8) then
             scale_l2g = view%scale_l2g_lookup(l2g_scale_type, view%lun_type(l))
             value = larr(l - beg + 1)
             if (value /= spval .and. scale_l2g /= spval) then
                sum_arr = sum_arr + value * scale_l2g * view%wtgcell_l(l)
                sumwt = sumwt + view%wtgcell_l(l)
             end if
          end if
       end do
       garr(g - g0 + 1) = weighted_average(sum_arr, sumwt)
    end do

  end subroutine l2g_range

  !-----------------------------------------------------------------------
  pure subroutine t2g_range(view, g0, g1, tarr, beg, garr)
    !
    ! !DESCRIPTION:
    ! For gridcells g0:g1, the average of each gridcell's topounits with
    ! nonzero weight whose value isn't spval (spval if there are none).
    ! Topounits have no active check and only the 'unity' scale type.
    !
    !$acc routine vector
    ! !ARGUMENTS:
    type(subgrid_ave_view), intent(in) :: view
    integer , intent(in)  :: g0, g1
    real(r8), intent(in)  :: tarr(*)          ! topounit values; tarr(1) is topounit beg
    integer , intent(in)  :: beg              ! index of tarr(1)
    real(r8), intent(out) :: garr(*)          ! garr(1) is gridcell g0
    !
    ! !LOCAL VARIABLES:
    integer  :: g, t
    real(r8) :: value, sum_arr, sumwt
    !-----------------------------------------------------------------------

    !$acc loop vector private(sum_arr, sumwt)
    do g = g0, g1
       sum_arr = 0._r8
       sumwt = 0._r8
       do t = view%topi(g), view%topf(g)
          value = tarr(t - beg + 1)
          if (view%wtgcell_t(t) /= 0._r8 .and. value /= spval) then
             sum_arr = sum_arr + value * view%wtgcell_t(t)
             sumwt = sumwt + view%wtgcell_t(t)
          end if
       end do
       garr(g - g0 + 1) = weighted_average(sum_arr, sumwt)
    end do

  end subroutine t2g_range

  !-----------------------------------------------------------------------
  pure real(r8) function weighted_average(sum_arr, sumwt)
    !
    ! !DESCRIPTION:
    ! sum_arr / sumwt, or spval if nothing contributed.
    !
    !$acc routine seq
    ! !ARGUMENTS:
    real(r8), intent(in) :: sum_arr, sumwt
    !-----------------------------------------------------------------------

    if (sumwt /= 0._r8) then
       weighted_average = sum_arr / sumwt
    else
       weighted_average = spval
    end if

  end function weighted_average

end module subgridAveViewMod
