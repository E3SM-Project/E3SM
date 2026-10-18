module fake_subgrid_mod

  !-----------------------------------------------------------------------
  ! Synthetic ELM subgrid for the history unit test: gridcells with 1-2
  ! topounits and 1-5 landunits, 1-4 columns per landunit and 1-4 patches per
  ! column, with proc bounds that don't start at 1. As in ELM (initGridCells),
  ! landunits are ordered by type across gridcells, and columns and patches
  ! follow their landunit, so a gridcell's points are not contiguous. Weights
  ! are products of random fractions (some zero), active flags and snow layers
  ! are random, and the c2l/l2g scale tables follow the shape of
  ! subgridAveMod's (factors for urban columns, spval for landunit types a
  ! scale type excludes). subgrid_view builds the subgrid_ave_view of it.
  !-----------------------------------------------------------------------

  use shr_kind_mod     , only : r8 => shr_kind_r8
  use subgridAveViewMod, only : subgrid_ave_view, build_child_map, &
       subgrid_ave_view_attach, spval
  use histComputeMod   , only : lev_grid, lev_topo, lev_lun, lev_col, lev_patch

  implicit none
  public

  integer , parameter :: nlevsno = 5
  integer , parameter :: n_c2l = 3, n_l2g = 6   ! unity/urbanf/urbans; unity/natveg/veg/ice/nonurb/lake
  integer , parameter :: n_ltype = 9, urb_min = 7

  integer :: begg, endg, begt, endt, begl, endl, begc, endc, begp, endp
  integer , allocatable, target :: grid_t(:), grid_l(:), grid_c(:), grid_p(:)
  integer , allocatable, target :: lun_c(:), lun_p(:), col_p(:), itype_l(:)
  real(r8), allocatable, target :: wt_t(:), wt_l(:), wt_c(:), wt_p(:)
  logical , allocatable, target :: act_t(:), act_l(:), act_c(:), act_p(:)
  integer , allocatable, target :: snl(:)
  real(r8), allocatable, target :: scale_c2l(:,:)          ! (column, c2l type)
  real(r8), allocatable, target :: scale_l2g_lookup(:,:)   ! (l2g type, landunit type)
  integer , allocatable, target :: grc_lun_ptr(:), grc_lun_idx(:)   ! landunits of a gridcell
  integer , allocatable, target :: topi(:), topf(:)                 ! topounits of a gridcell
  integer , allocatable, target :: coli(:), colf(:), pfti(:), pftf(:) ! columns, patches of a landunit
  real(r8), allocatable :: frac_l(:), frac_c(:), frac_p(:)   ! fractions of the parent
  logical :: extrasnowlayers = .false.

contains

  !-----------------------------------------------------------------------
  pure integer function level_beg(lev)
    integer, intent(in) :: lev
    select case (lev)
    case (lev_topo);  level_beg = begt
    case (lev_lun);   level_beg = begl
    case (lev_col);   level_beg = begc
    case (lev_patch); level_beg = begp
    case default;     level_beg = begg
    end select
  end function level_beg

  pure integer function level_end(lev)
    integer, intent(in) :: lev
    select case (lev)
    case (lev_topo);  level_end = endt
    case (lev_lun);   level_end = endl
    case (lev_col);   level_end = endc
    case (lev_patch); level_end = endp
    case default;     level_end = endg
    end select
  end function level_end

  real(r8) function urand()
    call random_number(urand)
  end function urand

  integer function rand_int(lo, hi)
    integer, intent(in) :: lo, hi
    rand_int = lo + min(hi - lo, int(urand() * (hi - lo + 1)))
  end function rand_int

  !-----------------------------------------------------------------------
  subroutine build_subgrid(ng, offset)
    integer, intent(in) :: ng      ! gridcells
    integer, intent(in) :: offset  ! proc bounds start at offset+1 (each level)
    integer, allocatable :: tg(:), lg(:), cg(:), cl(:), pg(:), pc(:), pl(:), lt(:)
    integer :: nt, nl, nc, np, g, i, ic, ip
    integer :: l, c, p, ltype, nl2, nc2, np2
    integer, allocatable :: lg2(:), lt2(:), cg2(:), cl2(:), pg2(:), pc2(:), pl2(:)

    allocate(tg(2*ng), lg(5*ng), lt(5*ng), cg(20*ng), cl(20*ng), pg(80*ng), pc(80*ng), pl(80*ng))
    nt = 0; nl = 0; nc = 0; np = 0
    do g = 1, ng
       do i = 1, rand_int(1, 2)
          nt = nt + 1; tg(nt) = g
       end do
       do i = 1, rand_int(1, 5)
          nl = nl + 1; lg(nl) = g; lt(nl) = rand_int(1, n_ltype)
          do ic = 1, rand_int(1, 4)
             nc = nc + 1; cg(nc) = g; cl(nc) = nl
             do ip = 1, rand_int(1, 4)
                np = np + 1; pg(np) = g; pc(np) = nc; pl(np) = nl
             end do
          end do
       end do
    end do

    ! reorder: landunits by type, then gridcell; columns and patches follow
    allocate(lg2(nl), lt2(nl), cg2(nc), cl2(nc), pg2(np), pc2(np), pl2(np))
    nl2 = 0; nc2 = 0; np2 = 0
    do ltype = 1, n_ltype
       do l = 1, nl
          if (lt(l) /= ltype) cycle
          nl2 = nl2 + 1; lg2(nl2) = lg(l); lt2(nl2) = lt(l)
          do c = 1, nc
             if (cl(c) /= l) cycle
             nc2 = nc2 + 1; cg2(nc2) = cg(c); cl2(nc2) = nl2
             do p = 1, np
                if (pc(p) /= c) cycle
                np2 = np2 + 1; pg2(np2) = pg(p); pc2(np2) = nc2; pl2(np2) = nl2
             end do
          end do
       end do
    end do
    lg(1:nl) = lg2; lt(1:nl) = lt2; cg(1:nc) = cg2; cl(1:nc) = cl2
    pg(1:np) = pg2; pc(1:np) = pc2; pl(1:np) = pl2

    begg = offset + 1; endg = offset + ng
    begt = offset + 1; endt = offset + nt
    begl = offset + 1; endl = offset + nl
    begc = offset + 1; endc = offset + nc
    begp = offset + 1; endp = offset + np
    allocate(grid_t(begt:endt), grid_l(begl:endl), grid_c(begc:endc), grid_p(begp:endp))
    allocate(lun_c(begc:endc), lun_p(begp:endp), col_p(begp:endp), itype_l(begl:endl))
    allocate(wt_t(begt:endt), wt_l(begl:endl), wt_c(begc:endc), wt_p(begp:endp))
    allocate(act_t(begt:endt), act_l(begl:endl), act_c(begc:endc), act_p(begp:endp))
    allocate(snl(begc:endc), frac_l(begl:endl), frac_c(begc:endc), frac_p(begp:endp))
    allocate(scale_c2l(begc:endc, n_c2l), scale_l2g_lookup(n_l2g, n_ltype))
    grid_t(:) = tg(1:nt) + offset
    grid_l(:) = lg(1:nl) + offset
    grid_c(:) = cg(1:nc) + offset
    grid_p(:) = pg(1:np) + offset
    lun_c(:)  = cl(1:nc) + offset
    lun_p(:)  = pl(1:np) + offset
    col_p(:)  = pc(1:np) + offset
    itype_l(:) = lt(1:nl)

    ! l2g: which landunit types each scale type keeps (init_scale_l2g_lookup)
    scale_l2g_lookup(1, :) = 1.0_r8
    do i = 2, n_l2g
       do ltype = 1, n_ltype
          scale_l2g_lookup(i, ltype) = merge(1.0_r8, spval, urand() < 0.5_r8)
       end do
    end do
    ! c2l: urban columns get canyon-geometry factors (init_scale_c2l)
    do c = begc, endc
       scale_c2l(c, 1) = 1.0_r8
       if (itype_l(lun_c(c)) >= urb_min) then
          scale_c2l(c, 2) = 3.0_r8 * (0.2_r8 + urand())
          scale_c2l(c, 3) = scale_c2l(c, 2) / (2._r8 * (0.2_r8 + urand()) + 1._r8)
       else
          scale_c2l(c, 2:3) = 1.0_r8
       end if
    end do

    call perturb_subgrid(1.0_r8)
  end subroutine build_subgrid

  !-----------------------------------------------------------------------
  subroutine perturb_subgrid(fraction)
    ! New weights / active flags / snow layers for a fraction of the points
    ! (1: all). Weights are products of fractions of the parent, as in ELM.
    real(r8), intent(in) :: fraction
    real(r8) :: s
    integer :: g, t, l, c, p

    do t = begt, endt
       if (urand() < fraction) wt_t(t) = merge(0._r8, urand(), urand() < 0.1_r8)
       if (urand() < fraction) act_t(t) = urand() < 0.9_r8
    end do
    do g = begg, endg      ! topounit weights of a gridcell sum to 1 (or 0)
       s = sum(wt_t, mask=grid_t == g)
       if (s > 0._r8) where (grid_t == g) wt_t = wt_t / s
    end do
    do l = begl, endl
       if (urand() < fraction) frac_l(l) = merge(0._r8, urand(), urand() < 0.1_r8)
       if (urand() < fraction) act_l(l) = urand() < 0.9_r8
    end do
    do g = begg, endg
       s = sum(frac_l, mask=grid_l == g)
       if (s > 0._r8) where (grid_l == g) frac_l = frac_l / s
    end do
    do c = begc, endc
       if (urand() < fraction) frac_c(c) = merge(0._r8, urand(), urand() < 0.1_r8)
       if (urand() < fraction) act_c(c) = urand() < 0.9_r8
       if (urand() < fraction) snl(c) = -rand_int(0, nlevsno)
    end do
    do l = begl, endl
       s = sum(frac_c, mask=lun_c == l)
       if (s > 0._r8) where (lun_c == l) frac_c = frac_c / s
    end do
    do p = begp, endp
       if (urand() < fraction) frac_p(p) = merge(0._r8, urand(), urand() < 0.1_r8)
       if (urand() < fraction) act_p(p) = urand() < 0.9_r8
    end do
    do c = begc, endc
       s = sum(frac_p, mask=col_p == c)
       if (s > 0._r8) where (col_p == c) frac_p = frac_p / s
    end do
    wt_l(:) = frac_l(:)
    do c = begc, endc
       wt_c(c) = wt_l(lun_c(c)) * frac_c(c)
    end do
    do p = begp, endp
       wt_p(p) = wt_c(col_p(p)) * frac_p(p)
    end do
    !$acc update device(wt_t, wt_l, wt_c, wt_p, act_t, act_l, act_c, act_p, snl) if_present
  end subroutine perturb_subgrid

  !-----------------------------------------------------------------------
  subroutine subgrid_view(view)
    ! The subgrid_ave_view of the synthetic subgrid, with ELM's index arrays
    ! (contiguous ranges and the landunits of each gridcell) built here, its
    ! arrays on the device and the view attached to them.
    type(subgrid_ave_view), intent(inout) :: view
    integer :: bad

    allocate(grc_lun_ptr(begg:endg+1), grc_lun_idx(begl:endl))
    call build_child_map(begg, endg, begl, endl, grid_l, grc_lun_ptr, grc_lun_idx, bad)
    if (bad /= 0) error stop 'landunit child map'
    allocate(topi(begg:endg), topf(begg:endg), coli(begl:endl), colf(begl:endl))
    allocate(pfti(begl:endl), pftf(begl:endl))
    call ranges(begt, begg, grid_t, topi, topf)
    call ranges(begc, begl, lun_c, coli, colf)
    call ranges(begp, begl, lun_p, pfti, pftf)

    view%begg = begg ; view%endg = endg ; view%begt = begt ; view%endt = endt
    view%begl = begl ; view%endl = endl ; view%begc = begc ; view%endc = endc
    view%begp = begp ; view%endp = endp
    view%grc_lun_ptr => grc_lun_ptr ; view%grc_lun_idx => grc_lun_idx
    view%topi => topi ; view%topf => topf
    view%coli => coli ; view%colf => colf ; view%pfti => pfti ; view%pftf => pftf
    view%wtgcell_t => wt_t ; view%wtgcell_l => wt_l ; view%wtgcell_c => wt_c ; view%wtgcell_p => wt_p
    view%active_t => act_t ; view%active_l => act_l ; view%active_c => act_c ; view%active_p => act_p
    view%col_of_pft => col_p ; view%lun_type => itype_l
    view%scale_c2l => scale_c2l ; view%scale_l2g_lookup => scale_l2g_lookup
    view%snl => snl

    !$acc enter data copyin(grc_lun_ptr, grc_lun_idx, topi, topf, coli, colf, pfti, pftf, &
    !$acc    wt_t, wt_l, wt_c, wt_p, act_t, act_l, act_c, act_p, col_p, itype_l, scale_c2l, &
    !$acc    scale_l2g_lookup, snl)
    call subgrid_ave_view_attach(view)

  contains

    subroutine ranges(first_child, first_parent, parent, first, last)
      ! first/last child of each parent; the children of a parent must be contiguous
      integer, intent(in)  :: first_child   ! index of parent(1)'s child
      integer, intent(in)  :: first_parent  ! lower bound of first/last
      integer, intent(in)  :: parent(:)     ! parent of each child
      integer, intent(out) :: first(first_parent:), last(first_parent:)
      integer :: i, child
      first(:) = huge(1)
      last(:)  = -huge(1)
      do i = 1, size(parent)
         child = first_child + i - 1
         first(parent(i)) = min(first(parent(i)), child)
         last(parent(i))  = max(last(parent(i)), child)
      end do
      do i = lbound(first,1), ubound(first,1)
         if (last(i) - first(i) + 1 /= count(parent == i)) error stop 'children not contiguous'
      end do
    end subroutine ranges

  end subroutine subgrid_view

end module fake_subgrid_mod
