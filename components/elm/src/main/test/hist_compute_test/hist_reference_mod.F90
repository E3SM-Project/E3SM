module hist_reference_mod

  !-----------------------------------------------------------------------
  ! Reference for the history unit test: the original ELM algorithm
  ! (histFileMod's hist_update_hbuf_field_1d/2d, hist_set_snow_field_2d and
  ! subgridAveMod's serial p2g/c2g/l2g/t2g), one field at a time, serially,
  ! on the test's synthetic subgrid, written without histComputeMod or
  ! subgridAveViewMod. It keeps the original rules for 'X' and 'M', and adds
  ! 'U'/'L' (extremes of the valid samples). One deliberate change, which
  ! histComputeMod makes too: the model array is indexed from its own lower
  ! bound for every flag (the original code did so only for 1d 'A' fields).
  !-----------------------------------------------------------------------

  use shr_kind_mod  , only : r8 => shr_kind_r8
  use histComputeMod, only : lev_grid, lev_topo, lev_lun, lev_col, lev_patch, &
       time_inst, time_mean, time_max, time_min, time_max_valid, time_min_valid, &
       no_snow_unset, no_snow_zero, hist_active_field, hist_src1d, hist_src2d
  use fake_subgrid_mod

  implicit none
  private

  public :: ref_update_field  ! accumulate one field
  public :: ref_normalize     ! hfields_normalize for one field
  public :: ref_p2g, ref_c2g, ref_l2g, ref_t2g  ! the serial averages to gridcells

contains

  !-----------------------------------------------------------------------
  subroutine ref_update_field(fld, src1, src2, hbuf, nacs)
    type(hist_active_field), intent(in)    :: fld
    type(hist_src1d)       , intent(in)    :: src1(:)
    type(hist_src2d)       , intent(in)    :: src2(:)
    real(r8)               , intent(inout) :: hbuf(:,:)   ! [buffer point (from 1), level]
    integer                , intent(inout) :: nacs(:,:)
    real(r8), allocatable :: field(:,:)     ! [point of the level (from beg), level]
    real(r8), allocatable :: field_gcell(:,:)
    logical , allocatable :: active(:)
    integer :: beg, end, j, k, nl, lb, lb2

    beg = level_beg(fld%source_level)
    end = level_end(fld%source_level)
    nl  = fld%num_levels
    allocate(field(beg:end, nl))
    if (fld%source_rank == 1) then
       lb = lbound(src1(fld%source_index)%ptr, 1)
       do k = beg, end
          field(k, 1) = src1(fld%source_index)%ptr(k - beg + lb)
       end do
    else
       ! levels 1:num2d of the array, except snow fields, whose levels
       ! hist_set_snow_field_2d counts from the array's lower bound
       lb = lbound(src2(fld%source_index)%ptr, 1)
       lb2 = 1
       if (fld%no_snow_behavior /= no_snow_unset) lb2 = lbound(src2(fld%source_index)%ptr, 2)
       do j = 1, nl
          do k = beg, end
             field(k, j) = src2(fld%source_index)%ptr(k - beg + lb, j - 1 + lb2)
          end do
       end do
    end if
    if (fld%no_snow_behavior /= no_snow_unset) call set_snow_field(fld, field)

    if (fld%average_to_gridcell .and. fld%source_level /= lev_grid) then
       allocate(field_gcell(begg:endg, nl))
       do j = 1, nl
          select case (fld%source_level)
          case (lev_patch)
             call ref_p2g(field(:, j), field_gcell(:, j), fld%c2l_scale, fld%l2g_scale)
          case (lev_col)
             call ref_c2g(field(:, j), field_gcell(:, j), fld%c2l_scale, fld%l2g_scale)
          case (lev_lun)
             call ref_l2g(field(:, j), field_gcell(:, j), fld%l2g_scale)
          case (lev_topo)
             call ref_t2g(field(:, j), field_gcell(:, j))
          end select
       end do
       call accumulate_gridcell(fld%time_op, field_gcell, hbuf, nacs)
    else
       allocate(active(beg:end))
       select case (fld%source_level)
       case (lev_topo);  active(:) = act_t(:)
       case (lev_lun);   active(:) = act_l(:)
       case (lev_col);   active(:) = act_c(:)
       case (lev_patch); active(:) = act_p(:)
       case default;     active(:) = .true.
       end select
       call accumulate_native(fld%time_op, field, active, hbuf, nacs)
    end if
  end subroutine ref_update_field

  !-----------------------------------------------------------------------
  subroutine accumulate_gridcell(time_op, field, hbuf, nacs)
    ! hist_update_hbuf_field_1d/2d, "map2gcell" branch (plus 'U'/'L')
    integer , intent(in)    :: time_op
    real(r8), intent(in)    :: field(:,:)
    real(r8), intent(inout) :: hbuf(:,:)
    integer , intent(inout) :: nacs(:,:)
    integer :: j, k

    do j = 1, size(field, 2)
       do k = 1, size(field, 1)
          select case (time_op)
          case (time_inst)
             if (field(k,j) /= spval) then
                hbuf(k,j) = field(k,j)
             else
                hbuf(k,j) = spval
             end if
             nacs(k,j) = 1
          case (time_mean)
             if (field(k,j) /= spval) then
                if (nacs(k,j) == 0) hbuf(k,j) = 0._r8
                hbuf(k,j) = hbuf(k,j) + field(k,j)
                nacs(k,j) = nacs(k,j) + 1
             else
                if (nacs(k,j) == 0) hbuf(k,j) = spval
             end if
          case (time_max)
             if (field(k,j) /= spval) then
                if (nacs(k,j) == 0) hbuf(k,j) = -1.e50_r8
                hbuf(k,j) = max(hbuf(k,j), field(k,j))
             else
                hbuf(k,j) = spval
             end if
             nacs(k,j) = 1
          case (time_min)
             if (field(k,j) /= spval) then
                if (nacs(k,j) == 0) hbuf(k,j) = +1.e50_r8
                hbuf(k,j) = min(hbuf(k,j), field(k,j))
             else
                hbuf(k,j) = spval
             end if
             nacs(k,j) = 1
          case (time_max_valid, time_min_valid)
             call accumulate_valid(time_op, field(k,j) /= spval, field(k,j), hbuf(k,j), nacs(k,j))
          end select
       end do
    end do
  end subroutine accumulate_gridcell

  !-----------------------------------------------------------------------
  subroutine accumulate_native(time_op, field, active, hbuf, nacs)
    ! hist_update_hbuf_field_1d/2d, native branch (plus 'U'/'L')
    integer , intent(in)    :: time_op
    real(r8), intent(in)    :: field(:,:)
    logical , intent(in)    :: active(:)
    real(r8), intent(inout) :: hbuf(:,:)
    integer , intent(inout) :: nacs(:,:)
    integer :: j, k

    do j = 1, size(field, 2)
       do k = 1, size(field, 1)
          select case (time_op)
          case (time_inst)
             if (active(k)) then
                if (field(k,j) /= spval) then
                   hbuf(k,j) = field(k,j)
                else
                   hbuf(k,j) = spval
                end if
             else
                hbuf(k,j) = spval
             end if
             nacs(k,j) = 1
          case (time_mean)
             if (active(k)) then
                if (field(k,j) /= spval) then
                   if (nacs(k,j) == 0) hbuf(k,j) = 0._r8
                   hbuf(k,j) = hbuf(k,j) + field(k,j)
                   nacs(k,j) = nacs(k,j) + 1
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
             else
                if (nacs(k,j) == 0) hbuf(k,j) = spval
             end if
          case (time_max)
             if (active(k)) then
                if (field(k,j) /= spval) then
                   if (nacs(k,j) == 0) hbuf(k,j) = -1.e50_r8
                   hbuf(k,j) = max(hbuf(k,j), field(k,j))
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
             else
                if (nacs(k,j) == 0) hbuf(k,j) = spval
             end if
             nacs(k,j) = 1
          case (time_min)
             if (active(k)) then
                if (field(k,j) /= spval) then
                   if (nacs(k,j) == 0) hbuf(k,j) = +1.e50_r8
                   hbuf(k,j) = min(hbuf(k,j), field(k,j))
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
             else
                if (nacs(k,j) == 0) hbuf(k,j) = spval
             end if
             nacs(k,j) = 1
          case (time_max_valid, time_min_valid)
             call accumulate_valid(time_op, active(k) .and. field(k,j) /= spval, field(k,j), &
                  hbuf(k,j), nacs(k,j))
          end select
       end do
    end do
  end subroutine accumulate_native

  !-----------------------------------------------------------------------
  subroutine accumulate_valid(time_op, valid, value, hbuf, nacs)
    ! 'U' / 'L': extremes of the valid samples
    integer , intent(in)    :: time_op
    logical , intent(in)    :: valid
    real(r8), intent(in)    :: value
    real(r8), intent(inout) :: hbuf
    integer , intent(inout) :: nacs

    if (valid) then
       if (time_op == time_max_valid) then
          if (nacs == 0) hbuf = -1.e50_r8
          hbuf = max(hbuf, value)
       else
          if (nacs == 0) hbuf = +1.e50_r8
          hbuf = min(hbuf, value)
       end if
       nacs = nacs + 1
    else
       if (nacs == 0) hbuf = spval
    end if
  end subroutine accumulate_valid

  !-----------------------------------------------------------------------
  subroutine ref_normalize(fld, hbuf, nacs)
    type(hist_active_field), intent(in)    :: fld
    real(r8)               , intent(inout) :: hbuf(:,:)
    integer                , intent(in)    :: nacs(:,:)
    integer :: j, k

    if (fld%time_op /= time_mean) return
    do j = 1, size(hbuf, 2)
       do k = 1, size(hbuf, 1)
          if (nacs(k,j) /= 0) hbuf(k,j) = hbuf(k,j) / float(nacs(k,j))
       end do
    end do
  end subroutine ref_normalize

  !-----------------------------------------------------------------------
  subroutine set_snow_field(fld, field)
    ! hist_set_snow_field_2d
    type(hist_active_field), intent(in)    :: fld
    real(r8)               , intent(inout) :: field(level_beg(fld%source_level):, :)
    real(r8), allocatable :: field_in(:,:)
    integer :: point, level, c, num_levels, num_snow_layers, num_nonexistent_layers
    real(r8) :: no_snow_val

    num_levels = size(field, 2)
    no_snow_val = spval
    if (fld%no_snow_behavior == no_snow_zero) no_snow_val = 0._r8
    allocate(field_in, source=field)
    do point = lbound(field, 1), ubound(field, 1)
       if (fld%source_level == lev_col) then
          c = point
       else
          c = col_p(point)
       end if
       num_snow_layers = abs(snl(c))
       num_nonexistent_layers = num_levels - num_snow_layers
       if (.not. extrasnowlayers) then
          do level = 1, num_nonexistent_layers
             field(point, level) = no_snow_val
          end do
          do level = (num_nonexistent_layers + 1), num_levels
             field(point, level) = field_in(point, level)
          end do
       else
          do level = num_levels, (num_levels-num_nonexistent_layers+1), -1
             field(point, level) = no_snow_val
          end do
          do level = (num_levels-num_nonexistent_layers), 1, -1
             field(point, level) = field_in(point, level+num_nonexistent_layers)
          end do
       end if
    end do
  end subroutine set_snow_field

  !-----------------------------------------------------------------------
  ! The original serial averages to gridcells: a scatter over the children
  ! in index order (p2c scale type is always 'unity').

  subroutine ref_p2g(parr, garr, ic2l, il2g)
    real(r8), intent(in)  :: parr(begp:)
    real(r8), intent(out) :: garr(begg:)
    integer , intent(in)  :: ic2l, il2g
    real(r8) :: sumwt(begg:endg), sc2l, sl2g
    integer :: p, c, l, g

    garr(:) = spval
    sumwt(:) = 0._r8
    do p = begp, endp
       if (act_p(p) .and. wt_p(p) /= 0._r8) then
          c = col_p(p)
          l = lun_p(p)
          sc2l = scale_c2l(c, ic2l)
          sl2g = scale_l2g_lookup(il2g, itype_l(l))
          if (parr(p) /= spval .and. sc2l /= spval .and. sl2g /= spval) then
             g = grid_p(p)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + parr(p) * sc2l * sl2g * wt_p(p)
             sumwt(g) = sumwt(g) + wt_p(p)
          end if
       end if
    end do
    call finish(garr, sumwt)
  end subroutine ref_p2g

  subroutine ref_c2g(carr, garr, ic2l, il2g)
    real(r8), intent(in)  :: carr(begc:)
    real(r8), intent(out) :: garr(begg:)
    integer , intent(in)  :: ic2l, il2g
    real(r8) :: sumwt(begg:endg), sc2l, sl2g
    integer :: c, l, g

    garr(:) = spval
    sumwt(:) = 0._r8
    do c = begc, endc
       if (act_c(c) .and. wt_c(c) /= 0._r8) then
          l = lun_c(c)
          sc2l = scale_c2l(c, ic2l)
          sl2g = scale_l2g_lookup(il2g, itype_l(l))
          if (carr(c) /= spval .and. sc2l /= spval .and. sl2g /= spval) then
             g = grid_c(c)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + carr(c) * sc2l * sl2g * wt_c(c)
             sumwt(g) = sumwt(g) + wt_c(c)
          end if
       end if
    end do
    call finish(garr, sumwt)
  end subroutine ref_c2g

  subroutine ref_l2g(larr, garr, il2g)
    real(r8), intent(in)  :: larr(begl:)
    real(r8), intent(out) :: garr(begg:)
    integer , intent(in)  :: il2g
    real(r8) :: sumwt(begg:endg), sl2g
    integer :: l, g

    garr(:) = spval
    sumwt(:) = 0._r8
    do l = begl, endl
       if (act_l(l) .and. wt_l(l) /= 0._r8) then
          sl2g = scale_l2g_lookup(il2g, itype_l(l))
          if (larr(l) /= spval .and. sl2g /= spval) then
             g = grid_l(l)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + larr(l) * sl2g * wt_l(l)
             sumwt(g) = sumwt(g) + wt_l(l)
          end if
       end if
    end do
    call finish(garr, sumwt)
  end subroutine ref_l2g

  subroutine ref_t2g(tarr, garr)
    real(r8), intent(in)  :: tarr(begt:)
    real(r8), intent(out) :: garr(begg:)
    real(r8) :: sumwt(begg:endg)
    integer :: t, g

    garr(:) = spval
    sumwt(:) = 0._r8
    do t = begt, endt
       if (wt_t(t) /= 0._r8 .and. tarr(t) /= spval) then
          g = grid_t(t)
          if (sumwt(g) == 0._r8) garr(g) = 0._r8
          garr(g) = garr(g) + tarr(t) * wt_t(t)
          sumwt(g) = sumwt(g) + wt_t(t)
       end if
    end do
    call finish(garr, sumwt)
  end subroutine ref_t2g

  subroutine finish(garr, sumwt)
    real(r8), intent(inout) :: garr(begg:)
    real(r8), intent(in)    :: sumwt(begg:)
    integer :: g

    do g = begg, endg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          error stop 'reference: sumwt is greater than 1.0'
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
  end subroutine finish

end module hist_reference_mod
