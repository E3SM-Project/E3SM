module hist_reference_mod

  !-----------------------------------------------------------------------
  ! Reference for the history unit test: the algorithm of histFileMod before
  ! histComputeMod (hist_update_hbuf_field_1d/2d, hist_set_snow_field_2d and
  ! subgridAveMod's p2g/c2g/l2g/t2g), one field at a time, serially, on the
  ! test's synthetic subgrid. Two deliberate changes, matching histComputeMod:
  ! - the model array is indexed from its own lower bound for every averaging
  !   flag (the old code did so only for 1d 'A' fields), and
  ! - 'X'/'M' use the 'A' rule: only valid samples count (nacs counts them) and
  !   a missing sample never overwrites the extreme found so far.
  !-----------------------------------------------------------------------

  use shr_kind_mod  , only : r8 => shr_kind_r8
  use histComputeMod, only : lev_grid, lev_topo, lev_lun, lev_col, lev_patch, &
       avg_inst, avg_mean, avg_max, avg_min, no_snow_unset, no_snow_zero, &
       hist_field_desc, hist_src1d, hist_src2d
  use fake_subgrid_mod

  implicit none
  private

  public :: ref_update_field  ! accumulate one field
  public :: ref_normalize     ! hfields_normalize for one field

contains

  !-----------------------------------------------------------------------
  subroutine ref_update_field(fd, src1, src2, hbuf, nacs)
    type(hist_field_desc), intent(in)    :: fd
    type(hist_src1d)     , intent(in)    :: src1(:)
    type(hist_src2d)     , intent(in)    :: src2(:)
    real(r8)             , intent(inout) :: hbuf(:,:)   ! [buffer point (from 1), level]
    integer              , intent(inout) :: nacs(:,:)
    real(r8), allocatable :: field(:,:)     ! [point of the level (from beg), level]
    real(r8), allocatable :: field_gcell(:,:)
    logical , allocatable :: active(:)
    integer :: beg, end, j, k, nl, lb, lb2

    beg = level_beg(fd%level)
    end = level_end(fd%level)
    nl  = fd%num2d
    allocate(field(beg:end, nl))
    if (fd%rank == 1) then
       lb = lbound(src1(fd%hpindex)%ptr, 1)
       do k = beg, end
          field(k, 1) = src1(fd%hpindex)%ptr(k - beg + lb)
       end do
    else
       ! levels 1:num2d of the array, except snow fields, whose levels
       ! hist_set_snow_field_2d counts from the array's lower bound
       lb = lbound(src2(fd%hpindex)%ptr, 1)
       lb2 = 1
       if (fd%no_snow /= no_snow_unset) lb2 = lbound(src2(fd%hpindex)%ptr, 2)
       do j = 1, nl
          do k = beg, end
             field(k, j) = src2(fd%hpindex)%ptr(k - beg + lb, j - 1 + lb2)
          end do
       end do
    end if
    if (fd%no_snow /= no_snow_unset) call set_snow_field(fd, field)

    if (fd%to_gridcell .and. fd%level /= lev_grid) then
       allocate(field_gcell(begg:endg, nl))
       do j = 1, nl
          select case (fd%level)
          case (lev_patch)
             call p2g(field(:, j), field_gcell(:, j), fd%c2l, fd%l2g)
          case (lev_col)
             call c2g(field(:, j), field_gcell(:, j), fd%c2l, fd%l2g)
          case (lev_lun)
             call l2g(field(:, j), field_gcell(:, j), fd%l2g)
          case (lev_topo)
             call t2g(field(:, j), field_gcell(:, j))
          end select
       end do
       allocate(active(begg:endg))
       active(:) = .true.
       call accumulate_all(fd%avg, field_gcell, active, hbuf, nacs)
    else
       allocate(active(beg:end))
       select case (fd%level)
       case (lev_topo);  active(:) = act_t(:)
       case (lev_lun);   active(:) = act_l(:)
       case (lev_col);   active(:) = act_c(:)
       case (lev_patch); active(:) = act_p(:)
       case default;     active(:) = .true.
       end select
       call accumulate_all(fd%avg, field, active, hbuf, nacs)
    end if
  end subroutine ref_update_field

  !-----------------------------------------------------------------------
  subroutine accumulate_all(avg, field, active, hbuf, nacs)
    integer , intent(in)    :: avg
    real(r8), intent(in)    :: field(:,:)
    logical , intent(in)    :: active(:)
    real(r8), intent(inout) :: hbuf(:,:)
    integer , intent(inout) :: nacs(:,:)
    integer :: j, k
    logical :: valid

    do j = 1, size(field, 2)
       do k = 1, size(field, 1)
          valid = active(k) .and. field(k,j) /= spval
          select case (avg)
          case (avg_inst)
             if (valid) then
                hbuf(k,j) = field(k,j)
             else
                hbuf(k,j) = spval
             end if
             nacs(k,j) = 1
          case (avg_mean)
             if (valid) then
                if (nacs(k,j) == 0) hbuf(k,j) = 0._r8
                hbuf(k,j) = hbuf(k,j) + field(k,j)
                nacs(k,j) = nacs(k,j) + 1
             else
                if (nacs(k,j) == 0) hbuf(k,j) = spval
             end if
          case (avg_max)
             if (valid) then
                if (nacs(k,j) == 0) hbuf(k,j) = -1.e50_r8
                hbuf(k,j) = max(hbuf(k,j), field(k,j))
                nacs(k,j) = nacs(k,j) + 1
             else
                if (nacs(k,j) == 0) hbuf(k,j) = spval
             end if
          case (avg_min)
             if (valid) then
                if (nacs(k,j) == 0) hbuf(k,j) = +1.e50_r8
                hbuf(k,j) = min(hbuf(k,j), field(k,j))
                nacs(k,j) = nacs(k,j) + 1
             else
                if (nacs(k,j) == 0) hbuf(k,j) = spval
             end if
          end select
       end do
    end do
  end subroutine accumulate_all

  !-----------------------------------------------------------------------
  subroutine ref_normalize(fd, hbuf, nacs)
    type(hist_field_desc), intent(in)    :: fd
    real(r8)             , intent(inout) :: hbuf(:,:)
    integer              , intent(in)    :: nacs(:,:)
    integer :: j, k

    if (fd%avg /= avg_mean) return
    do j = 1, size(hbuf, 2)
       do k = 1, size(hbuf, 1)
          if (nacs(k,j) /= 0) hbuf(k,j) = hbuf(k,j) / float(nacs(k,j))
       end do
    end do
  end subroutine ref_normalize

  !-----------------------------------------------------------------------
  subroutine set_snow_field(fd, field)
    ! hist_set_snow_field_2d
    type(hist_field_desc), intent(in)    :: fd
    real(r8)             , intent(inout) :: field(level_beg(fd%level):, :)
    real(r8), allocatable :: field_in(:,:)
    integer :: point, level, c, num_levels, num_snow_layers, num_nonexistent_layers
    real(r8) :: no_snow_val

    num_levels = size(field, 2)
    no_snow_val = spval
    if (fd%no_snow == no_snow_zero) no_snow_val = 0._r8
    allocate(field_in, source=field)
    do point = lbound(field, 1), ubound(field, 1)
       if (fd%level == lev_col) then
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
  ! subgridAveMod's averages to gridcells (the scale tables are those built
  ! by create_scale_c2l / build_scale_l2g; p2c is always 'unity')

  subroutine p2g(parr, garr, ic2l, il2g)
    real(r8), intent(in)  :: parr(begp:)
    real(r8), intent(out) :: garr(begg:)
    integer , intent(in)  :: ic2l, il2g
    real(r8) :: scale_p2c(begp:endp), sumwt(begg:endg), sc2l(begc:endc), sl2g(begl:endl)
    integer :: p, c, l, g

    call create_scales(ic2l, il2g, sc2l, sl2g)
    scale_p2c(:) = 1.0_r8
    garr(:) = spval
    sumwt(:) = 0._r8
    do p = begp, endp
       if (act_p(p) .and. wt_p(p) /= 0._r8) then
          c = col_p(p)
          l = lun_p(p)
          if (parr(p) /= spval .and. sc2l(c) /= spval .and. sl2g(l) /= spval) then
             g = grid_p(p)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + parr(p) * scale_p2c(p) * sc2l(c) * sl2g(l) * wt_p(p)
             sumwt(g) = sumwt(g) + wt_p(p)
          end if
       end if
    end do
    call finish(garr, sumwt)
  end subroutine p2g

  subroutine c2g(carr, garr, ic2l, il2g)
    real(r8), intent(in)  :: carr(begc:)
    real(r8), intent(out) :: garr(begg:)
    integer , intent(in)  :: ic2l, il2g
    real(r8) :: sumwt(begg:endg), sc2l(begc:endc), sl2g(begl:endl)
    integer :: c, l, g

    call create_scales(ic2l, il2g, sc2l, sl2g)
    garr(:) = spval
    sumwt(:) = 0._r8
    do c = begc, endc
       if (act_c(c) .and. wt_c(c) /= 0._r8) then
          l = lun_c(c)
          if (carr(c) /= spval .and. sc2l(c) /= spval .and. sl2g(l) /= spval) then
             g = grid_c(c)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + carr(c) * sc2l(c) * sl2g(l) * wt_c(c)
             sumwt(g) = sumwt(g) + wt_c(c)
          end if
       end if
    end do
    call finish(garr, sumwt)
  end subroutine c2g

  subroutine l2g(larr, garr, il2g)
    real(r8), intent(in)  :: larr(begl:)
    real(r8), intent(out) :: garr(begg:)
    integer , intent(in)  :: il2g
    real(r8) :: sumwt(begg:endg), sc2l(begc:endc), sl2g(begl:endl)
    integer :: l, g

    call create_scales(1, il2g, sc2l, sl2g)
    garr(:) = spval
    sumwt(:) = 0._r8
    do l = begl, endl
       if (act_l(l) .and. wt_l(l) /= 0._r8) then
          if (larr(l) /= spval .and. sl2g(l) /= spval) then
             g = grid_l(l)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + larr(l) * sl2g(l) * wt_l(l)
             sumwt(g) = sumwt(g) + wt_l(l)
          end if
       end if
    end do
    call finish(garr, sumwt)
  end subroutine l2g

  subroutine t2g(tarr, garr)
    real(r8), intent(in)  :: tarr(begt:)
    real(r8), intent(out) :: garr(begg:)
    real(r8) :: scale_t2g(begt:endt), sumwt(begg:endg)
    integer :: t, g

    scale_t2g(:) = 1.0_r8
    garr(:) = spval
    sumwt(:) = 0._r8
    do t = begt, endt
       if (wt_t(t) /= 0._r8) then
          if (tarr(t) /= spval .and. scale_t2g(t) /= spval) then
             g = grid_t(t)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + tarr(t) * scale_t2g(t) * wt_t(t)
             sumwt(g) = sumwt(g) + wt_t(t)
          end if
       end if
    end do
    call finish(garr, sumwt)
  end subroutine t2g

  subroutine create_scales(ic2l, il2g, sc2l, sl2g)
    ! create_scale_c2l / build_scale_l2g: the scale factors of all columns and
    ! landunits for one scale type, rebuilt for every average as before
    integer , intent(in)  :: ic2l, il2g
    real(r8), intent(out) :: sc2l(begc:), sl2g(begl:)
    integer :: c, l
    do c = begc, endc
       sc2l(c) = scale_c2l(c, ic2l)
    end do
    do l = begl, endl
       sl2g(l) = scale_l2g(l, il2g)
    end do
  end subroutine create_scales

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
